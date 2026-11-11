"""Validated triangle geometry and intrinsic transport in physical millimeters."""

from dataclasses import dataclass
import numpy as np


@dataclass
class Geometry:
    points: np.ndarray
    faces: np.ndarray
    frames: np.ndarray
    xy: np.ndarray
    areas: np.ndarray
    neighbors: np.ndarray
    neighbor_edges: np.ndarray
    transport: np.ndarray


def vertex_normals(points, faces):
    tri = points[faces]
    normal = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    out = np.zeros_like(points, dtype=float)
    for corner in range(3):
        np.add.at(out, faces[:, corner], normal)
    lengths = np.linalg.norm(out, axis=1)
    if np.any(lengths < 1e-14):
        raise ValueError("Cannot compute normals at isolated or degenerate vertices")
    return out / lengths[:, None]


def build_geometry(points, faces):
    points = np.ascontiguousarray(points, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Surface points must be finite N x 3 world coordinates")
    if faces.ndim != 2 or faces.shape[1] != 3 or not len(faces):
        raise ValueError("A nonempty triangular surface is required")
    if faces.min() < 0 or faces.max() >= len(points):
        raise ValueError("Invalid triangle vertex index")
    tri = points[faces]
    e1 = tri[:, 1]-tri[:, 0]
    cross = np.cross(e1, tri[:, 2]-tri[:, 0])
    magnitude = np.linalg.norm(cross, axis=1)
    if np.any(magnitude <= 1e-12):
        raise ValueError("Degenerate triangle; repair the mesh before tracking")
    normal = cross/magnitude[:, None]
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(normal, e1)
    frames = np.ascontiguousarray(np.stack((e1, e2, normal), axis=2))
    xy = np.ascontiguousarray(np.einsum('fvi,fij->fvj', tri-tri[:, :1], frames[:, :, :2]))
    neighbors = np.full((len(faces), 3), -1, dtype=np.int64)
    neighbor_edges = np.full_like(neighbors, -1)
    transport = np.zeros((len(faces), 3, 2, 2))
    a = faces.ravel()
    b = np.roll(faces, -1, axis=1).ravel()
    low, high = np.minimum(a, b), np.maximum(a, b)
    order = np.lexsort((high, low))
    breaks = np.r_[True, (low[order[1:]] != low[order[:-1]]) | (high[order[1:]] != high[order[:-1]]), True]
    starts = np.flatnonzero(breaks)
    multiplicity = np.diff(starts)
    if np.any(multiplicity > 2):
        raise ValueError("Nonmanifold edge has more than two incident faces")
    paired = starts[:-1][multiplicity == 2]
    left, right = order[paired], order[paired+1]
    if np.any(a[left] != b[right]) or np.any(b[left] != a[right]):
        raise ValueError("Inconsistent face winding across a shared edge")
    f, e, g, q = left//3, left % 3, right//3, right % 3
    axis = points[b[left]]-points[a[left]]
    axis /= np.linalg.norm(axis, axis=1)[:, None]
    cosine = np.clip(np.einsum('ni,ni->n', normal[f], normal[g]), -1, 1)
    sine = np.einsum('ni,ni->n', axis, np.cross(normal[f], normal[g]))
    angle = np.arctan2(sine, cosine)
    cosine, sine = np.cos(angle), np.sin(angle)
    skew = np.zeros((len(axis), 3, 3))
    skew[:, 0, 1], skew[:, 0, 2] = -axis[:, 2], axis[:, 1]
    skew[:, 1, 0], skew[:, 1, 2] = axis[:, 2], -axis[:, 0]
    skew[:, 2, 0], skew[:, 2, 1] = -axis[:, 1], axis[:, 0]
    rotation = cosine[:, None, None]*np.eye(3)+(1-cosine)[:, None, None]*axis[:, :, None]*axis[:, None, :]+sine[:, None, None]*skew
    local = np.einsum('nai,nab,nbj->nij', frames[g, :, :2], rotation, frames[f, :, :2])
    if not np.allclose(np.einsum('nji,njk->nik', local, local), np.eye(2), atol=1e-9):
        raise ValueError("Nonorthogonal intrinsic edge transport")
    neighbors[f, e], neighbors[g, q] = g, f
    neighbor_edges[f, e], neighbor_edges[g, q] = q, e
    transport[f, e], transport[g, q] = local, local.transpose(0, 2, 1)
    return Geometry(points, faces, frames, xy, magnitude/2, neighbors, neighbor_edges, np.ascontiguousarray(transport))


def restrict_surface(points, faces, point_data, mask_name="roi"):
    if mask_name not in point_data:
        mask = np.ones(len(points), dtype=bool)
    else:
        mask = np.asarray(point_data[mask_name]).ravel() > 0
    kept_faces = np.all(mask[faces], axis=1)
    if not kept_faces.any():
        raise ValueError("ROI contains no complete triangles")
    # Retain all ROI vertex IDs, matching the legacy mask convention.
    ids = np.flatnonzero(mask)
    mapping = np.full(len(points), -1, dtype=np.int64)
    mapping[ids] = np.arange(len(ids))
    return points[ids], mapping[faces[kept_faces]], {name: values[ids] for name, values in point_data.items()}, ids, np.flatnonzero(kept_faces)
