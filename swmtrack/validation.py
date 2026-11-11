"""Independent binary format checks and physical-coordinate validation."""

from pathlib import Path
import struct
import re
import hashlib
import numpy as np
import nibabel as nib
from scipy.spatial import cKDTree


def streamlines_digest(streamlines):
    digest = hashlib.sha256()
    for points in streamlines:
        digest.update(struct.pack('<q', len(points)))
        digest.update(np.asarray(points, dtype='<f4').tobytes())
    return digest.hexdigest()


def read_tck_independent(path):
    with open(path, 'rb') as stream:
        header = b''
        while not header.endswith(b'END\n'):
            line = stream.readline()
            if not line or len(header) > 100000:
                raise ValueError('Invalid TCK header')
            header += line
        offset = int(re.search(rb'file: \. (\d+)', header).group(1))
        if b'datatype: Float32LE' not in header:
            raise ValueError('Expected little-endian float32 TCK')
        stream.seek(offset)
        rows = np.frombuffer(stream.read(), dtype='<f4').reshape(-1, 3)
    if not len(rows) or not np.isposinf(rows[-1]).all():
        raise ValueError('TCK infinity EOF marker is missing')
    out, first = [], 0
    for k in np.flatnonzero(np.isnan(rows).all(axis=1)):
        out.append(rows[first:k].copy()); first = k+1
    if first != len(rows)-1:
        raise ValueError('Invalid TCK streamline separators')
    return out


def read_trk_independent(path):
    with open(path, 'rb') as stream:
        header = stream.read(1000)
        if header[:6] != b'TRACK\0' or struct.unpack_from('<i', header, 996)[0] != 1000:
            raise ValueError('Invalid TRK header')
        if struct.unpack_from('<h', header, 36)[0] or struct.unpack_from('<h', header, 238)[0]:
            raise ValueError('This validator expects no TRK scalar/property payload')
        spacing = np.array(struct.unpack_from('<3f', header, 12))
        affine = np.array(struct.unpack_from('<16f', header, 440)).reshape(4, 4)
        count = struct.unpack_from('<i', header, 988)[0]
        out = []
        for _ in range(count):
            n = struct.unpack('<i', stream.read(4))[0]
            points = np.frombuffer(stream.read(12*n), dtype='<f4').reshape(n, 3)
            out.append(nib.affines.apply_affine(affine, points/spacing-.5))
        if stream.read():
            raise ValueError('Unexpected trailing TRK data')
    return out


def point_face_error(geometry, points, face_ids):
    q = points-geometry.points[geometry.faces[face_ids, 0]]
    local = np.einsum('ni,nij->nj', q, geometry.frames[face_ids])
    b2 = local[:, 1]/geometry.xy[face_ids, 2, 1]
    b1 = (local[:, 0]-geometry.xy[face_ids, 2, 0]*b2)/geometry.xy[face_ids, 1, 0]
    bary = np.column_stack((1-b1-b2, b1, b2))
    return np.abs(local[:, 2]), np.min(bary, axis=1)


def point_triangle_distance(geometry, points, face_ids):
    """Euclidean distance to each closed triangle, in physical millimeters.

    Barycentric errors grow inversely with triangle altitude. A fixed
    barycentric cutoff can therefore reject a correctly rounded float32
    endpoint on a thin triangle. Use the existing physical tolerance instead.
    """
    distances, bary = point_face_error(geometry, points, face_ids)
    outside = bary < 0
    if np.any(outside):
        triangles = geometry.points[geometry.faces[face_ids[outside]]]
        query = np.asarray(points)[outside]
        edge_distance = np.full(len(query), np.inf)
        for corner in range(3):
            start = triangles[:, corner]
            edge = triangles[:, (corner+1) % 3]-start
            fraction = np.clip(np.einsum('ni,ni->n', query-start, edge)/
                               np.einsum('ni,ni->n', edge, edge), 0, 1)
            closest = start+fraction[:, None]*edge
            edge_distance = np.minimum(edge_distance, np.linalg.norm(query-closest, axis=1))
        distances[outside] = edge_distance
    return distances


def surface_point_distances(geometry, points, tolerance=1e-4, candidates=32):
    """Locate surface points with a bounded fallback for long/thin triangles.

    Nearest triangle centers are only an acceleration heuristic. Every triangle
    lies inside its vertex bounding sphere, so the fallback includes every
    triangle that could be within the requested physical tolerance.
    """
    points=np.asarray(points)
    centers=geometry.points[geometry.faces].mean(axis=1)
    radii=np.linalg.norm(geometry.points[geometry.faces]-centers[:,None,:],axis=2).max(axis=1)
    tree=cKDTree(centers)
    _,near=tree.query(points,k=min(candidates,len(geometry.faces)))
    if near.ndim==1:near=near[:,None]
    best=np.full(len(points),np.inf);best_face=np.full(len(points),-1,np.int64)
    for j in range(near.shape[1]):
        distances=point_triangle_distance(geometry,points,near[:,j])
        improved=distances<best;best[improved]=distances[improved];best_face[improved]=near[improved,j]
    fallback=np.flatnonzero(best>=tolerance)
    for i in fallback:
        ids=np.asarray(tree.query_ball_point(points[i],float(radii.max()+tolerance)),np.int64)
        eligible=np.linalg.norm(centers[ids]-points[i],axis=1)<=radii[ids]+tolerance
        ids=ids[eligible]
        if not len(ids):continue
        distances=point_triangle_distance(geometry,np.repeat(points[i:i+1],len(ids),axis=0),ids)
        nearest=int(np.argmin(distances))
        if distances[nearest]<best[i]:best[i]=distances[nearest];best_face[i]=ids[nearest]
    return best,best_face,len(fallback)


def validate_exports(prefix, streamlines, reference, prepared=None, diagnostics=None):
    report = {'streamlines':len(streamlines), 'points':int(sum(len(s) for s in streamlines)),
              'streamline_digest':streamlines_digest(streamlines)}
    for suffix, reader in [('.tck', read_tck_independent), ('.trk', read_trk_independent)]:
        loaded = reader(str(prefix)+suffix)
        assert len(loaded) == len(streamlines)
        error = max((float(np.max(np.abs(x-y))) for x,y in zip(loaded,streamlines)), default=0.)
        assert error < 3e-5, (suffix,error)
        nib_loaded = nib.streamlines.load(str(prefix)+suffix)
        assert len(nib_loaded.streamlines) == len(streamlines)
        report[suffix[1:]+'_independent_roundtrip_max_abs_mm'] = error
    points = np.concatenate(streamlines) if streamlines else np.empty((0,3))
    assert np.isfinite(points).all()
    voxel = nib.affines.apply_affine(np.linalg.inv(reference.affine), points)
    assert np.all(voxel >= -1e-4) and np.all(voxel <= np.array(reference.shape[:3])-1+1e-4)
    report['all_points_finite_and_inside_fod'] = True
    if prepared is not None and len(points):
        g = prepared.geometry
        sample = points[np.linspace(0,len(points)-1,min(len(points),20000),dtype=int)]
        best,nearest_faces,fallback_count=surface_point_distances(g,sample)
        assert np.isfinite(best).all() and best.max() < 1e-4, float(best.max())
        plane,_=point_face_error(g,sample,nearest_faces)
        report['surface_check'] = {'sampled_points':len(sample),'max_plane_error_mm':float(plane.max()),
            'max_triangle_distance_mm':float(best.max()),'bounded_search_fallback_points':fallback_count,
            'distance_metric':'Euclidean distance to the closed triangle; center search with bounding-sphere fallback'}
        if diagnostics is not None and len(streamlines):
            ids = diagnostics['accepted_trial_ids']
            endpoint_faces = diagnostics['endpoint_faces'][ids].reshape(-1)
            endpoints = np.asarray([[s[0],s[-1]] for s in streamlines]).reshape(-1,3)
            assert np.all(endpoint_faces >= 0) and np.all(endpoint_faces < len(g.faces))
            distances = point_triangle_distance(g,endpoints,endpoint_faces)
            assert np.max(distances) < 1e-4, float(np.max(distances))
            report['endpoint_face_check_passed'] = True
            report['endpoint_face_max_distance_mm'] = float(np.max(distances))
            report['endpoint_face_tolerance_mm'] = 1e-4
            flags = diagnostics['stops'][ids]
            anatomical = np.isin(flags, [1, 2])
            policy = str(diagnostics.get('endpoint_policy', 'unknown'))
            if policy == 'within_gyrus':
                if prepared.gyrus is None:
                    raise ValueError('within_gyrus validation requires the unified gyrus mask')
                target_labels = prepared.gyrus
            else:
                target_labels = prepared.targets
            actual_labels = target_labels[endpoint_faces].reshape(-1, 2)
            assert np.array_equal(actual_labels[anatomical], flags[anatomical])
            if 'endpoint_labels' in diagnostics:
                assert np.array_equal(actual_labels, diagnostics['endpoint_labels'][ids])
            report['anatomical_endpoint_labels_match'] = True
            report['all_endpoints_in_distinct_gyri'] = bool(np.all(np.sort(actual_labels, axis=1) == [1, 2]))
            if policy == 'between_gyri':
                assert report['all_endpoints_in_distinct_gyri']
            if policy == 'within_gyrus':
                report['all_endpoints_in_unified_gyrus'] = bool(np.all(actual_labels == 1))
                assert report['all_endpoints_in_unified_gyrus']
    return report
