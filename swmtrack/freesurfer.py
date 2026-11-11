"""Import native FreeSurfer geometry without changing vertex identities.

Surface RAS (tkregister RAS), scanner RAS, and FOD world RAS are distinct.
An identity scanner-to-world mapping must be explicitly requested; anatomical
registration cannot be established from matching image headers alone.
"""

from pathlib import Path
import json
import numpy as np
import nibabel as nib
from nibabel.freesurfer.io import read_geometry, read_annot, read_morph_data, read_label
from .geometry import vertex_normals
from .io import sha256, write_surface, read_surface


def affine_matrix(value):
    matrix = np.asarray(value, dtype=float)
    if (matrix.shape != (4, 4) or not np.isfinite(matrix).all()
            or not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-12)
            or abs(np.linalg.det(matrix[:3, :3])) < 1e-12):
        raise ValueError("Transform must be a finite, invertible 4 x 4 affine")
    return matrix


def surface_ras_to_scanner(reference):
    """FreeSurfer voxel-to-scanner-RAS @ inverse(voxel-to-tkregister-RAS)."""
    if not hasattr(reference.header, "get_vox2ras_tkr"):
        raise ValueError("A FreeSurfer MGH/MGZ reference is required")
    return affine_matrix(reference.affine @ np.linalg.inv(reference.header.get_vox2ras_tkr()))


def import_freesurfer(subject_dir, reference_fod, output_dir, *,
                      scanner_to_world=None, assume_scanner_aligned=False,
                      hemispheres=("lh", "rh"), annotation="aparc",
                      gyrus_a=None, gyrus_b=None):
    """Write white/pial VTKs, native atlas IDs, morphometry, and provenance.

    ``scanner_to_world`` acts on physical RAS millimeters, not voxel indices.
    Nonlinear transforms are deliberately not approximated by an affine here.
    Full-gyrus masks are optional and are not legacy gyral skeleton masks.
    """
    if (scanner_to_world is None) == (not assume_scanner_aligned):
        raise ValueError("Supply scanner_to_world OR explicitly set assume_scanner_aligned=True")
    if not hemispheres or any(h not in ("lh", "rh") for h in hemispheres):
        raise ValueError("Select lh and/or rh")
    if not annotation or Path(annotation).name != annotation:
        raise ValueError("annotation must be a FreeSurfer annotation name")
    root, output = Path(subject_dir).resolve(), Path(output_dir).resolve()
    reference = nib.load(root / "mri/orig.mgz")
    fod = nib.load(reference_fod)
    if len(fod.shape) != 4:
        raise ValueError("FOD must be a 4D SH coefficient image")
    source_to_scanner = surface_ras_to_scanner(reference)
    scanner_to_world = affine_matrix(np.eye(4) if scanner_to_world is None else scanner_to_world)
    source_to_world = scanner_to_world @ source_to_scanner
    inverse_fod = np.linalg.inv(fod.affine)
    metadata = {
        "coordinate_system": "RAS+", "units": "mm", "subject_dir": str(root),
        "reference_fod": str(Path(reference_fod).resolve()),
        "fod_shape": list(fod.shape), "fod_affine": fod.affine.tolist(),
        "surface_ras_to_scanner_ras": source_to_scanner.tolist(),
        "scanner_ras_to_fod_world": scanner_to_world.tolist(),
        "surface_ras_to_fod_world": source_to_world.tolist(),
        "registration": "scanner alignment explicitly assumed; requires anatomical QC"
                        if assume_scanner_aligned else "supplied physical RAS affine; requires anatomical QC",
        "annotation": annotation, "endpoint_mask_definition": "full annotation regions, not gyral skeletons",
        "source_hashes": {"mri/orig.mgz": sha256(root / "mri/orig.mgz")},
        "hemispheres": {},
    }
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / "surface_ras_to_fod_world.txt", source_to_world, fmt="%.12g")
    np.savetxt(output / "scanner_ras_to_fod_world.txt", scanner_to_world, fmt="%.12g")
    for hemi in hemispheres:
        white, faces, surface_meta = read_geometry(str(root / "surf" / f"{hemi}.white"), read_metadata=True)
        pial, pial_faces = read_geometry(str(root / "surf" / f"{hemi}.pial"))
        if white.shape != pial.shape or not np.array_equal(faces, pial_faces):
            raise ValueError(f"{hemi}: white and pial topology differs")
        if (not np.isfinite(white).all() or not np.isfinite(pial).all()
                or faces.ndim != 2 or faces.shape[1] != 3
                or faces.min() < 0 or faces.max() >= len(white)):
            raise ValueError(f"{hemi}: invalid surface geometry")
        if "cras" in surface_meta and not np.allclose(surface_meta["cras"], source_to_scanner[:3, 3], atol=1e-3):
            raise ValueError(f"{hemi}: surface volume center differs from orig.mgz; check reference provenance")
        labels, table, names = read_annot(str(root / "label" / f"{hemi}.{annotation}.annot"))
        if len(labels) != len(white) or np.any(labels < -1) or np.any(labels >= len(names)):
            raise ValueError(f"{hemi}: annotation does not match surface vertices")
        names = [name.decode("utf-8") for name in names]
        atlas = np.zeros(len(labels), dtype=np.int32)
        assigned = labels >= 0
        atlas[assigned] = table[labels[assigned], -1]
        arrays = {"source_vertex_id": np.arange(len(white), dtype=np.int32),
                  annotation: atlas, "annotation_index": labels.astype(np.int32)}
        sources = [f"surf/{hemi}.white", f"surf/{hemi}.pial", f"label/{hemi}.{annotation}.annot"]
        for field in ("curv", "sulc", "thickness"):
            relative = f"surf/{hemi}.{field}"
            if (root / relative).exists():
                data = read_morph_data(str(root / relative))
                if len(data) != len(white) or not np.isfinite(data).all():
                    raise ValueError(f"{hemi}: invalid {field} values")
                arrays[field] = data
                sources.append(relative)
        cortex_file = f"label/{hemi}.cortex.label"
        if (root / cortex_file).exists():
            indices = read_label(str(root / cortex_file))
            if np.any(indices < 0) or np.any(indices >= len(white)):
                raise ValueError(f"{hemi}: invalid cortex vertex IDs")
            cortex = np.zeros(len(white), dtype=np.uint8); cortex[indices] = 1
            arrays["cortex"] = cortex
            sources.append(cortex_file)
        for key, name in (("gyrus_a", gyrus_a), ("gyrus_b", gyrus_b)):
            if name is not None:
                if name not in names:
                    raise ValueError(f"{hemi}: annotation has no region named {name}")
                arrays[key] = (labels == names.index(name)).astype(np.uint8)
        for relative in sources:
            metadata["source_hashes"][relative] = sha256(root / relative)
        # Reflection changes handedness; keep outward normals and winding consistent.
        if np.linalg.det(source_to_world[:3, :3]) < 0:
            faces = faces[:, [0, 2, 1]]
        hemisphere = {"vertices": len(white), "faces": len(faces), "surfaces": {},
                      "labels": [{"index": i, "name": name, "id": int(table[i, -1]),
                                  "rgba": [*map(int, table[i, :3]), 255-int(table[i, 3])],
                                  "vertices": int(np.count_nonzero(labels == i))}
                                 for i, name in enumerate(names)]}
        for kind, original in (("white", white), ("pial", pial)):
            points = nib.affines.apply_affine(source_to_world, original)
            data = {**arrays, "normals": vertex_normals(points, faces)}
            path = output / f"{hemi}.{kind}.world.vtk"
            write_surface(path, points, faces, data)
            rp, rf, rd = read_surface(path)
            if not np.array_equal(rf, faces) or any(not np.array_equal(rd[k], a) for k, a in data.items()):
                raise RuntimeError("VTK round trip changed topology or point arrays")
            error = float(np.max(np.abs(rp-points)))
            if error > 1e-6:
                raise RuntimeError("VTK round trip changed world coordinates")
            voxel = nib.affines.apply_affine(inverse_fod, points)
            inside = np.all((voxel >= -1e-5) & (voxel <= np.asarray(fod.shape[:3])-1+1e-5), axis=1)
            hemisphere["surfaces"][kind] = {
                "path": str(path), "sha256": sha256(path), "vtk_roundtrip_max_error_mm": error,
                "in_fod_domain_fraction": float(inside.mean()),
                "world_bounds_mm": [points.min(0).tolist(), points.max(0).tolist()],
            }
        metadata["hemispheres"][hemi] = hemisphere
    destination = output / "freesurfer_import.json"
    destination.write_text(json.dumps(metadata, indent=2)+"\n")
    return metadata
