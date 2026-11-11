"""Preparation, versioned caches, area-based seeding, and bounded tracking batches."""

from dataclasses import dataclass, asdict
from pathlib import Path
from collections import Counter
import json
import time
import numpy as np
import nibabel as nib
import numba
from . import ALGORITHM_VERSION, __version__
from .geometry import Geometry, build_geometry, restrict_surface, vertex_normals
from .harmonics import degree_from_count, project_coefficients
from .io import read_surface, write_surface, sample_sh, sha256
from .kernel import track_batch, STOP_NAMES


@dataclass
class Prepared:
    geometry: Geometry
    coefficients: np.ndarray
    seed_mask: np.ndarray
    targets: np.ndarray
    metadata: dict
    gyrus: np.ndarray | None = None


@dataclass
class TrackingConfig:
    seed: int = 42
    threads: int = 4
    seed_mode: str = "area"
    n_seeds: int | None = None
    seeds_per_mm2: float = 20.0
    samples_per_face: int = 6
    step_mm: float = 0.5
    max_turn_deg_per_mm: float = 14.0
    angular_resolution_deg: float = 2.0
    min_fod: float = 0.01
    min_length_mm: float = 8.0
    max_length_mm: float = 150.0
    endpoint_policy: str = "between_gyri"
    stop_at_targets: bool | None = None
    max_points_per_half: int = 4096
    batch_size: int = 128

    def validate(self):
        numeric = (self.step_mm, self.min_length_mm, self.max_length_mm, self.max_turn_deg_per_mm,
                   self.angular_resolution_deg, self.min_fod, self.seeds_per_mm2)
        if not np.isfinite(numeric).all():
            raise ValueError("Tracking parameters must be finite")
        if self.endpoint_policy not in ("between_gyri", "within_gyrus", "unrestricted"):
            raise ValueError("Endpoint policy must be between_gyri, within_gyrus, or unrestricted")
        if self.seed_mode not in ("area", "centroid"):
            raise ValueError("Seed mode must be area or centroid")
        if self.step_mm <= 0 or self.min_length_mm < 0 or self.max_length_mm <= self.min_length_mm:
            raise ValueError("Invalid physical step or tract length limits")
        if self.max_turn_deg_per_mm < 0 or self.max_turn_deg_per_mm*self.step_mm > 180:
            raise ValueError("Physical turning limit must give an angle in [0, 180]")
        if not 0 < self.angular_resolution_deg <= 180 or self.min_fod < 0:
            raise ValueError("Invalid angular resolution or FOD threshold")
        if self.threads < 1 or self.batch_size < 1 or self.samples_per_face < 1 or self.max_points_per_half < 2:
            raise ValueError("Invalid worker, batch, seed, or point limit")
        if self.seed < 0 or self.seed >= 2**63 or self.seeds_per_mm2 <= 0:
            raise ValueError("Invalid RNG seed or seed density")
        if self.n_seeds is not None and self.n_seeds < 1:
            raise ValueError("n_seeds must be positive")
        constrained = self.endpoint_policy in ("between_gyri", "within_gyrus")
        stop = constrained if self.stop_at_targets is None else self.stop_at_targets
        return stop


def unified_gyrus_faces(values, faces, label=None):
    """Select ONE anatomical label; never silently union several atlas gyri."""
    values = np.asarray(values).ravel()
    if not np.isfinite(values).all():
        raise ValueError("Gyrus labels must be finite")
    if label is None:
        labels = np.unique(values[values > 0])
        if len(labels) > 1:
            raise ValueError("Gyrus array contains multiple labels; select one with --gyrus-label")
        mask = values > 0
    else:
        if label <= 0:
            raise ValueError("Gyrus label must be positive; zero is background")
        mask = values == label
    return np.ascontiguousarray(mask[faces].any(axis=1), dtype=np.int64)


def endpoint_targets(prepared, policy):
    if policy == "within_gyrus":
        if prepared.gyrus is None or not np.any(prepared.gyrus):
            raise ValueError("within_gyrus requires one nonempty gyrus label; prepare with a 'gyrus' point array or --gyrus-array")
        if prepared.gyrus.shape != (len(prepared.geometry.faces),) or not np.isin(prepared.gyrus, [0, 1]).all():
            raise ValueError("Unified gyrus faces must be a binary face mask")
        return np.ascontiguousarray(prepared.gyrus, dtype=np.int64)
    return prepared.targets


def prepare(surface_path, fod_path, shrink_mm=0.5, interpolation="vertex", sh_to_world=None, roi_output=None,
            gyrus_array="gyrus", gyrus_label=None):
    start = time.perf_counter()
    if interpolation not in ("vertex", "centroid") or shrink_mm < 0:
        raise ValueError("Interpolation must be vertex or centroid; shrink distance must be nonnegative")
    points, faces, data = read_surface(surface_path)
    normals = np.asarray(data["normals"] if "normals" in data else vertex_normals(points, faces), dtype=float)
    normal_lengths = np.linalg.norm(normals, axis=1)
    if not np.isfinite(normal_lengths).all() or np.any(normal_lengths < 1e-12):
        raise ValueError("Surface normals are invalid")
    points = points - shrink_mm*normals/normal_lengths[:, None]
    points, faces, data, ids, face_ids = restrict_surface(points, faces, data)
    if gyrus_array in data:
        gyrus = unified_gyrus_faces(data[gyrus_array], faces, gyrus_label)
    elif gyrus_label is not None:
        raise ValueError(f"Gyrus point array {gyrus_array!r} is missing")
    else:
        gyrus = None
    loaded = time.perf_counter()
    geometry = build_geometry(points, faces)
    geometric = time.perf_counter()
    image = nib.load(fod_path)
    if len(image.shape) != 4:
        raise ValueError("FOD must be four-dimensional")
    degree = degree_from_count(image.shape[3])
    orientation = np.eye(3) if sh_to_world is None else np.asarray(sh_to_world, dtype=float)
    if orientation.shape != (3, 3) or not np.allclose(orientation.T@orientation, np.eye(3), atol=1e-7):
        raise ValueError("SH-to-world direction transform must be an orthogonal 3 x 3 matrix")
    local_frames = np.einsum('ij,fjk->fik', orientation.T, geometry.frames)
    if np.linalg.det(orientation) < 0:
        local_frames *= -1  # Even-degree SH are invariant to full inversion.
    count = 3 if interpolation == "vertex" else 1
    coefficients = np.empty((len(faces), count, degree+1), dtype=np.float32)
    if interpolation == "vertex":
        sampled = sample_sh(image, points)
    else:
        sampled = sample_sh(image, points[faces].mean(axis=1))
    sampled_time = time.perf_counter()
    for first in range(0, len(faces), 4096):
        last = min(first+4096, len(faces))
        for corner in range(count):
            values = sampled[faces[first:last, corner]] if count == 3 else sampled[first:last]
            coefficients[first:last, corner] = project_coefficients(values, local_frames[first:last], degree)
    projected = time.perf_counter()
    seed_vertices = np.asarray(data.get("seed", np.ones(len(points)))).ravel() > 0
    seed_mask = seed_vertices[faces].any(axis=1)
    a = np.asarray(data.get("gyrus_a", np.zeros(len(points)))).ravel() > 0
    b = np.asarray(data.get("gyrus_b", np.zeros(len(points)))).ravel() > 0
    targets = a[faces].any(axis=1).astype(np.int64) + 2*b[faces].any(axis=1).astype(np.int64)
    stat = Path(fod_path).stat()
    metadata = {"algorithm_version": ALGORITHM_VERSION, "package_version": __version__,
        "coordinate_system": "RAS+", "units": "mm", "sh_degree": degree,
        "sh_basis": "even-l_m-ascending_CS_negative-sine_positive-cosine",
        "projection": "exact_linear_polar_marginal_fourier_then_clip_in_2D",
        "sh_to_world": orientation.tolist(), "interpolation": interpolation, "shrink_mm": shrink_mm,
        "reference_fod": str(Path(fod_path).resolve()), "reference_affine": image.affine.tolist(),
        "reference_shape": list(image.shape), "fod_sha256": sha256(fod_path),
        "fod_stat": [stat.st_size, stat.st_mtime_ns], "surface_sha256": sha256(surface_path),
        "roi_vertices": len(points), "roi_triangles": len(faces), "seed_triangles": int(seed_mask.sum()),
        "seed_area_mm2": float(geometry.areas[seed_mask].sum()), "overlapping_target_faces": int((targets == 3).sum()),
        "target_overlap_policy": "stop_as_ambiguous_and_reject_between_gyri",
        "unified_gyrus_array": gyrus_array if gyrus is not None else None,
        "unified_gyrus_label": gyrus_label,
        "unified_gyrus_faces": int(gyrus.sum()) if gyrus is not None else 0,
        "timing_seconds": {"surface_load": loaded-start, "geometry": geometric-loaded,
                           "sample_sh": sampled_time-geometric, "project_coefficients": projected-sampled_time}}
    if roi_output:
        data["source_vertex_id"] = data.get("source_vertex_id", ids)
        write_surface(roi_output, points, faces, data)
    return Prepared(geometry, np.ascontiguousarray(coefficients), np.asarray(seed_mask), np.ascontiguousarray(targets), metadata, gyrus)


def save_prepared(path, prepared):
    arrays = {"geometry_"+name: value for name, value in vars(prepared.geometry).items()}
    arrays.update(coefficients=prepared.coefficients, seed_mask=prepared.seed_mask, targets=prepared.targets,
                  metadata=np.array(json.dumps(prepared.metadata)))
    if prepared.gyrus is not None:
        arrays["gyrus"] = prepared.gyrus
    np.savez(path, **arrays)


def load_prepared(path):
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata"]))
        if metadata["algorithm_version"] != ALGORITHM_VERSION:
            raise ValueError("Cache algorithm version mismatch; prepare again")
        geometry = Geometry(**{name: data["geometry_"+name].copy() for name in Geometry.__dataclass_fields__})
        coefficients, seeds, targets = [data[name].copy() for name in ("coefficients", "seed_mask", "targets")]
        gyrus = data["gyrus"].copy() if "gyrus" in data else None
    return Prepared(geometry, coefficients, seeds, targets, metadata, gyrus)


def generate_seeds(prepared, config):
    rng = np.random.default_rng(config.seed)
    eligible = np.flatnonzero(prepared.seed_mask)
    if not len(eligible):
        raise ValueError("No seed triangles")
    if config.seed_mode == "centroid":
        faces = np.tile(eligible, config.samples_per_face)
        if config.n_seeds is not None:
            faces = rng.choice(faces, size=config.n_seeds, replace=config.n_seeds > len(faces))
        barycentric = np.full((len(faces), 3), 1/3)
    else:
        areas = prepared.geometry.areas[eligible]
        count = config.n_seeds if config.n_seeds is not None else max(1, int(np.ceil(areas.sum()*config.seeds_per_mm2)))
        faces = rng.choice(eligible, size=count, p=areas/areas.sum())
        u, v = rng.random((2, count)); root = np.sqrt(u)
        barycentric = np.column_stack((1-root, root*(1-v), root*v))
    points = np.einsum('ni,nij->nj', barycentric, prepared.geometry.points[prepared.geometry.faces[faces]])
    return np.ascontiguousarray(faces), np.ascontiguousarray(points)


def track(prepared, config=None):
    config = config or TrackingConfig()
    stop_at_targets = config.validate()
    if config.endpoint_policy == "between_gyri" and (not np.any(prepared.targets & 1) or not np.any(prepared.targets & 2)):
        raise ValueError("between_gyri requires gyrus_a and gyrus_b point labels")
    targets = endpoint_targets(prepared, config.endpoint_policy)
    g = prepared.geometry
    seed_faces, seed_points = generate_seeds(prepared, config)
    stops = np.zeros((len(seed_faces), 2), dtype=np.int64)
    endpoint_faces = np.full_like(stops, -1)
    endpoint_labels = np.zeros_like(stops)
    lengths = np.zeros(len(seed_faces))
    accepted = np.zeros(len(seed_faces), dtype=bool)
    operations = np.zeros((len(seed_faces), 2), dtype=np.int64)
    streamlines, accepted_ids = [], []
    original_threads = numba.get_num_threads()
    if config.threads > numba.config.NUMBA_NUM_THREADS:
        raise ValueError("Requested worker count exceeds the Numba thread pool")
    numba.set_num_threads(config.threads)
    start = time.perf_counter(); batch_seconds = []
    try:
        for first in range(0, len(seed_faces), config.batch_size):
            last = min(first+config.batch_size, len(seed_faces))
            tick = time.perf_counter()
            paths, counts, flags, ends, ll, ops = track_batch(g.points, g.faces, g.frames, g.xy, g.neighbors,
                g.transport, prepared.coefficients, targets, seed_faces[first:last], seed_points[first:last],
                np.arange(first, last, dtype=np.int64), config.seed, config.step_mm,
                np.deg2rad(config.max_turn_deg_per_mm*config.step_mm), np.deg2rad(config.angular_resolution_deg),
                config.min_fod, config.max_length_mm, config.max_points_per_half, stop_at_targets)
            batch_seconds.append(time.perf_counter()-tick)
            stops[first:last], endpoint_faces[first:last], lengths[first:last], operations[first:last] = flags, ends, ll, ops
            labels = np.zeros_like(ends)
            valid_endpoints = ends >= 0
            labels[valid_endpoints] = targets[ends[valid_endpoints]]
            endpoint_labels[first:last] = labels
            good = (counts >= 2) & (ll >= config.min_length_mm) & (ll <= config.max_length_mm+1e-8)
            good &= ~np.any((flags == 6) | (flags == 8), axis=1)
            if config.endpoint_policy == "between_gyri":
                good &= ((labels[:, 0] == 1) & (labels[:, 1] == 2)) | ((labels[:, 0] == 2) & (labels[:, 1] == 1))
            elif config.endpoint_policy == "within_gyrus":
                good &= (labels[:, 0] == 1) & (labels[:, 1] == 1)
            accepted[first:last] = good
            for local in np.flatnonzero(good):
                streamlines.append(np.array(paths[local, :counts[local]], dtype=np.float32))
                accepted_ids.append(first+int(local))
    finally:
        numba.set_num_threads(original_threads)
    elapsed = time.perf_counter()-start
    stop_names = dict(STOP_NAMES)
    if config.endpoint_policy == "within_gyrus":
        stop_names[1] = "gyrus"
    report = {"algorithm_version": ALGORITHM_VERSION, "package_version": __version__, "config": asdict(config), "stop_at_targets": bool(stop_at_targets),
        "attempts": len(seed_faces), "accepted": len(streamlines), "acceptance_rate": len(streamlines)/len(seed_faces),
        "tracking_seconds_including_JIT_if_cold": elapsed, "first_batch_seconds": batch_seconds[0],
        "remaining_batches_seconds": float(sum(batch_seconds[1:])),
        "stop_counts": {stop_names[int(k)]: int(v) for k, v in Counter(stops.ravel()).items()},
        "integration_samples": int(operations[:, 0].sum()), "triangle_crossings": int(operations[:, 1].sum()),
        "stored_points": int(sum(len(s) for s in streamlines)),
        "length_mm": {str(q): float(np.percentile(lengths[accepted], q)) for q in (0, 50, 95, 100)} if accepted.any() else {}}
    diagnostics = {"stops": stops, "endpoint_faces": endpoint_faces, "lengths_mm": lengths,
        "endpoint_labels": endpoint_labels,
        "accepted": accepted, "accepted_trial_ids": np.asarray(accepted_ids, dtype=np.int64),
        "seed_faces": seed_faces, "seed_points": seed_points, "operations": operations,
        "endpoint_policy": np.array(config.endpoint_policy)}
    return streamlines, report, diagnostics
