"""Standard world-space surface and tractogram I/O plus explicit legacy import."""

from pathlib import Path
import json
import hashlib
import meshio
import nibabel as nib
import numpy as np
from .geometry import vertex_normals


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_legacy_obj(path):
    """MNI/BIC ASCII polygon OBJ (not Wavefront OBJ), with 0-based indices."""
    with open(path) as stream:
        if stream.read(1) != "P":
            raise ValueError("Expected MNI/BIC polygon OBJ starting with P")
        data = np.fromstring(stream.read(), sep=" ")
    n = int(data[5]); offset = 6
    points = data[offset:offset+3*n].reshape(n, 3).copy(); offset += 3*n
    normals = data[offset:offset+3*n].reshape(n, 3).copy(); offset += 3*n
    count, color_type = int(data[offset]), int(data[offset+1]); offset += 2
    if color_type not in (0, 1, 2):
        raise ValueError("Unsupported MNI polygon color encoding")
    offset += 4 * (1 if color_type == 0 else count if color_type == 1 else n)
    ends = data[offset:offset+count].astype(np.int64); offset += count
    if not np.array_equal(ends, np.arange(1, count+1)*3):
        raise ValueError("Only triangular MNI polygon surfaces are supported")
    faces = data[offset:offset+3*count].reshape(count, 3).astype(np.int64)
    if offset+3*count != len(data):
        raise ValueError("Unexpected trailing MNI OBJ data")
    return points, faces, normals


def read_surface(path):
    mesh = meshio.read(path)
    if any(block.type != "triangle" for block in mesh.cells):
        raise ValueError("Tracking requires triangle cells only")
    faces = np.concatenate([block.data for block in mesh.cells])
    return np.asarray(mesh.points, dtype=float), np.asarray(faces, dtype=np.int64), dict(mesh.point_data)


def write_surface(path, points, faces, point_data=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    data = {name: np.asarray(value, dtype=np.uint8) if np.asarray(value).dtype == bool else np.asarray(value) for name, value in (point_data or {}).items()}
    mesh = meshio.Mesh(np.asarray(points), [("triangle", np.asarray(faces))], point_data=data)
    if path.suffix.lower() == ".vtk":
        # VTK 5.1 emits vtktypeint32/vtktypeuint8, which older ParaView
        # LegacyVTKReader versions reject. VTK 4.2 uses portable type names.
        meshio.write(path, mesh, file_format="vtk42", binary=True)
    else:
        meshio.write(path, mesh)


def convert_legacy(mesh_path, fod_path, output, masks=None, input_space="image-mm", mesh_to_world=None):
    points, faces, normals = read_legacy_obj(mesh_path)
    image = nib.load(fod_path)
    if mesh_to_world is not None:
        transform = np.asarray(mesh_to_world, dtype=float)
    elif input_space == "world":
        transform = np.eye(4)
    elif input_space == "image-mm":
        transform = image.affine @ np.diag([*(1/np.asarray(image.header.get_zooms()[:3])), 1])
    elif input_space == "voxel":
        transform = image.affine
    else:
        raise ValueError("Input space must be world, image-mm, or voxel; supply an explicit mesh-to-world matrix for other spaces")
    if transform.shape != (4, 4) or not np.isfinite(transform).all():
        raise ValueError("mesh-to-world must be a finite 4 x 4 transform")
    points = nib.affines.apply_affine(transform, points)
    normals = normals @ np.linalg.inv(transform[:3, :3])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    if np.linalg.det(transform[:3, :3]) < 0:
        faces = faces[:, [0, 2, 1]]
    data = {"normals": normals, "source_vertex_id": np.arange(len(points), dtype=np.int64)}
    for name, path in (masks or {}).items():
        values = np.fromfile(path, dtype="<i4")
        if len(values) != len(points):
            raise ValueError(f"Mask {name} has {len(values)} entries, expected {len(points)}")
        # Preserve unified/atlas gyrus IDs so preparation can select ONE label.
        data[name] = values.copy() if name == "gyrus" else (values > 0).astype(np.uint8)
    write_surface(output, points, faces, data)
    metadata = {"coordinate_system": "RAS+", "units": "mm", "input_space": input_space,
                "mesh_to_world": transform.tolist(), "source_mesh": str(Path(mesh_path).resolve()),
                "source_mesh_sha256": sha256(mesh_path), "reference_fod": str(Path(fod_path).resolve())}
    Path(str(output)+".json").write_text(json.dumps(metadata, indent=2)+"\n")
    return metadata


def sample_sh(image, world_points):
    """Trilinear sampling with full affine, bounds, datatype, and NIfTI scaling.

    Uncompressed inputs use a read-only memory map; only sampled coefficients
    are copied. Compressed input is decompressed once, not once per coefficient.
    """
    voxel = nib.affines.apply_affine(np.linalg.inv(image.affine), world_points)
    if len(image.shape) != 4 or any(s < 2 for s in image.shape[:3]):
        raise ValueError("FOD must be a 4D coefficient image with spatial dimensions >= 2")
    upper = np.asarray(image.shape[:3])-1
    # NIfTI-1 stores affine rows in float32; allow their boundary roundoff.
    boundary_tolerance = 1e-5
    if np.any(voxel < -boundary_tolerance) or np.any(voxel > upper+boundary_tolerance):
        bad = np.flatnonzero(np.any((voxel < -boundary_tolerance) | (voxel > upper+boundary_tolerance), axis=1))
        raise ValueError(f"{len(bad)} surface samples are outside the FOD domain")
    voxel = np.clip(voxel, 0, upper)
    lower = np.minimum(np.floor(voxel).astype(np.int64), upper.astype(int)-1)
    fraction = voxel-lower
    proxy = image.dataobj
    raw = proxy.get_unscaled() if hasattr(proxy, "get_unscaled") else np.asanyarray(proxy)
    slope = float(getattr(proxy, "slope", 1))
    intercept = float(getattr(proxy, "inter", 0))
    out = np.zeros((len(voxel), image.shape[3]), dtype=np.float64)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                weight = np.prod(np.where(np.array([dx, dy, dz]), fraction, 1-fraction), axis=1)
                values = raw[lower[:, 0]+dx, lower[:, 1]+dy, lower[:, 2]+dz, :]
                out += weight[:, None]*(np.asarray(values, dtype=float)*slope+intercept)
    if not np.isfinite(out).all():
        raise ValueError("Nonfinite interpolated SH coefficients")
    return out


def write_streamlines_vtk(path,streamlines,cell_data=None):
    """Legacy binary POLYDATA with bounded connectivity batches and cell scalars."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields={}
    for name,values in (cell_data or {}).items():
        values=np.asarray(values)
        if not name or any(c.isspace() for c in name) or values.shape!=(len(streamlines),):
            raise ValueError('VTK cell scalars need a simple name and one value per track')
        if not np.isfinite(values).all():raise ValueError('VTK cell scalars must be finite')
        if values.dtype.kind in 'iu':
            if len(values) and (values.min()<np.iinfo(np.int32).min or values.max()>np.iinfo(np.int32).max):
                raise ValueError('VTK integer scalar exceeds int32 range')
            fields[name]=(values.astype('>i4'),'int')
        else:fields[name]=(values.astype('>f8'),'double')
    n=sum(len(p) for p in streamlines)
    if n+len(streamlines)>np.iinfo(np.int32).max:raise ValueError('VTK connectivity exceeds int32 capacity')
    with path.open('wb') as stream:
        stream.write(b'# vtk DataFile Version 4.2\nSWM physical RAS mm streamlines\nBINARY\nDATASET POLYDATA\n')
        stream.write(f'POINTS {n} float\n'.encode())
        # Avoid one small filesystem write per trajectory on shared storage.
        for first in range(0,len(streamlines),4096):
            np.concatenate(streamlines[first:first+4096]).astype('>f4').tofile(stream)
        stream.write(f'\nLINES {len(streamlines)} {n+len(streamlines)}\n'.encode())
        offset=0
        for first in range(0,len(streamlines),4096):
            block=[]
            for p in streamlines[first:first+4096]:
                block.append(np.r_[len(p),np.arange(offset,offset+len(p))]);offset+=len(p)
            np.concatenate(block).astype('>i4').tofile(stream)
        stream.write(b'\n')
        if fields:
            stream.write(f'CELL_DATA {len(streamlines)}\n'.encode())
            for name,(values,kind) in fields.items():
                stream.write(f'SCALARS {name} {kind} 1\nLOOKUP_TABLE default\n'.encode())
                values.tofile(stream);stream.write(b'\n')


def write_tractograms(prefix, streamlines, reference, vtk=False,*,vtk_cell_data=None):
    prefix = Path(prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    tractogram = nib.streamlines.Tractogram(streamlines, affine_to_rasmm=np.eye(4))
    tck = Path(str(prefix)+".tck"); trk = Path(str(prefix)+".trk")
    nib.streamlines.save(nib.streamlines.TckFile(tractogram), tck)
    header = nib.streamlines.TrkFile.create_empty_header()
    header["dimensions"] = reference.shape[:3]
    header["voxel_sizes"] = nib.affines.voxel_sizes(reference.affine)
    header["voxel_to_rasmm"] = reference.affine
    header["voxel_order"] = "".join(nib.aff2axcodes(reference.affine)).encode("ascii")
    nib.streamlines.save(nib.streamlines.TrkFile(tractogram, header=header), trk)
    result={"tck": str(tck), "trk": str(trk)}
    if vtk:
        path=Path(str(prefix)+'.vtk')
        write_streamlines_vtk(path,streamlines,vtk_cell_data)
        result['vtk']=str(path)
    return result
