"""Independent curvature estimators for an unchanged, indexed white surface.

Both estimators work on a smoothed copy of the vertices and return one feature
per original vertex ID. Neither method moves the tracking/segmentation mesh.
"""

import heapq

import numpy as np
from numba import njit
from scipy import ndimage
from scipy.spatial import cKDTree

from .geometry import vertex_normals


@njit(cache=True)
def _quadric_at_vertices(points, normals, area, indptr, indices, lengths,
                         radius, maximum_neighbors, sources):
    n = len(points)
    low = np.empty(len(sources), dtype=np.float64)
    high = np.empty(len(sources), dtype=np.float64)
    used = np.zeros(n, dtype=np.int32)
    distance = np.empty(n, dtype=np.float64)
    neighborhood = np.empty(maximum_neighbors, dtype=np.int64)
    path_lengths = np.empty(maximum_neighbors, dtype=np.float64)
    expanded_count = 0
    maximum_effective_radius = radius
    for output_index in range(len(sources)):
        source = sources[output_index]
        current_stamp = output_index + 1
        used[source] = current_stamp
        distance[source] = 0.
        queue = [(0., source)]
        count = 0
        while queue:
            d, v = heapq.heappop(queue)
            if d > radius and count >= 9:
                break
            if d > distance[v] + 1e-12:
                continue
            if count == maximum_neighbors:
                raise ValueError('Quadric neighborhood exceeded maximum_neighbors')
            neighborhood[count] = v
            path_lengths[count] = d
            count += 1
            for edge in range(indptr[v], indptr[v + 1]):
                neighbor = np.int64(indices[edge])
                candidate = d + lengths[edge]
                if used[neighbor] != current_stamp or candidate < distance[neighbor] - 1e-12:
                    used[neighbor] = current_stamp
                    distance[neighbor] = candidate
                    heapq.heappush(queue, (candidate, neighbor))
        if count < 9:
            raise ValueError('Quadric fit needs at least nine vertices in its connected component')
        effective_radius = max(radius, path_lengths[count - 1])
        if effective_radius > radius:
            expanded_count += 1
            maximum_effective_radius = max(maximum_effective_radius, effective_radius)
        scale = effective_radius / 2
        normal = normals[source]
        axis = np.zeros(3)
        axis[0 if abs(normal[0]) < .9 else 1] = 1.
        tangent1 = np.cross(axis, normal)
        tangent1 /= np.sqrt(np.sum(tangent1 * tangent1))
        tangent2 = np.cross(normal, tangent1)
        ata = np.zeros((6, 6))
        atb = np.zeros(6)
        for item in range(count):
            v = neighborhood[item]
            delta = points[v] - points[source]
            x = np.dot(delta, tangent1)
            y = np.dot(delta, tangent2)
            z = np.dot(delta, normal)
            basis = np.array([1., x, y, .5*x*x, x*y, .5*y*y])
            weight = area[v] * np.exp(-.5 * (path_lengths[item] / scale) ** 2)
            for a in range(6):
                atb[a] += weight * basis[a] * z
                for b in range(6):
                    ata[a, b] += weight * basis[a] * basis[b]
        ridge = np.trace(ata) * 1e-10 / 6
        for a in range(6):
            ata[a, a] += ridge
        coefficients = np.linalg.solve(ata, atb)
        slope_x, slope_y = coefficients[1], coefficients[2]
        a, b, c = coefficients[3], coefficients[4], coefficients[5]
        q = 1 + slope_x*slope_x + slope_y*slope_y
        mean = ((1+slope_y*slope_y)*a - 2*slope_x*slope_y*b +
                (1+slope_x*slope_x)*c) / (2*q**1.5)
        gaussian = (a*c-b*b) / (q*q)
        spread = np.sqrt(max(mean*mean-gaussian, 0.))
        low[output_index] = mean - spread
        high[output_index] = mean + spread
    return low, high, expanded_count, maximum_effective_radius


def quadric_principal_curvatures(points, faces, radius_mm=3., maximum_neighbors=512,
                                 source_ids=None):
    """Area-weighted Monge fits over connected geodesic neighborhoods.

    The local fit includes offset and linear terms so small normal errors do
    not masquerade as curvature. A physical graph radius prevents a fit from
    crossing a narrow sulcus merely because two banks are close in space.
    """
    from .open_preprocessing import mesh_operators, weighted_graph
    if radius_mm <= 0 or maximum_neighbors < 9:
        raise ValueError('Quadric radius and maximum_neighbors must be positive')
    _, area, _ = mesh_operators(points, faces)
    graph, _, _ = weighted_graph(points, faces)
    normals = vertex_normals(points, faces)
    if source_ids is None:
        source_ids = np.arange(len(points), dtype=np.int64)
    else:
        source_ids = np.asarray(source_ids, dtype=np.int64)
        if source_ids.ndim != 1 or np.any(source_ids < 0) or np.any(source_ids >= len(points)):
            raise ValueError('Invalid source_ids for quadric curvature')
    low, high, expanded_count, maximum_effective_radius = _quadric_at_vertices(np.asarray(points, dtype=np.float64),
                                     normals, area, graph.indptr, graph.indices,
                                     graph.data, float(radius_mm), int(maximum_neighbors), source_ids)
    if not np.isfinite(low).all() or not np.isfinite(high).all():
        raise ValueError('Nonfinite quadric curvature')
    return low, high, {'radius_mm': float(radius_mm),
                       'maximum_neighbors': int(maximum_neighbors),
                       'expanded_neighborhood_count': int(expanded_count),
                       'maximum_effective_radius_mm': float(maximum_effective_radius)}


@njit(cache=True)
def _closest_on_segment(p, a, b):
    ab = b - a
    length2 = np.dot(ab, ab)
    t = max(0., min(1., np.dot(p - a, ab) / length2))
    q = a + t * ab
    return q, np.dot(p - q, p - q)


@njit(cache=True)
def _closest_on_triangle(p, a, b, c, unit_normal):
    ab, ac = b - a, c - a
    pa = p - a
    signed_plane = np.dot(pa, unit_normal)
    projected = p - signed_plane * unit_normal
    vector = projected - a
    aa, bb, cc = np.dot(ab, ab), np.dot(ab, ac), np.dot(ac, ac)
    denominator = aa * cc - bb * bb
    u = (cc * np.dot(vector, ab) - bb * np.dot(vector, ac)) / denominator
    v = (aa * np.dot(vector, ac) - bb * np.dot(vector, ab)) / denominator
    if u >= 0. and v >= 0. and u + v <= 1.:
        return projected, signed_plane * signed_plane
    q, best = _closest_on_segment(p, a, b)
    for x, y in ((b, c), (c, a)):
        other, distance2 = _closest_on_segment(p, x, y)
        if distance2 < best:
            q, best = other, distance2
    return q, best


@njit(cache=True)
def _signed_distances(query, candidates, triangles, triangle_normals):
    result = np.empty(len(query), dtype=np.float32)
    for i in range(len(query)):
        nearest2 = np.inf
        signed = 0.
        for j in range(candidates.shape[1]):
            face = candidates[i, j]
            q, distance2 = _closest_on_triangle(query[i], triangles[face, 0],
                                                 triangles[face, 1], triangles[face, 2],
                                                 triangle_normals[face])
            if distance2 < nearest2:
                nearest2 = distance2
                signed = np.dot(query[i] - q, triangle_normals[face])
        result[i] = np.sqrt(nearest2) * (1. if signed >= 0. else -1.)
    return result


def _surface_signed_distance_grid(points, faces, spacing_mm, band_mm,
                                  candidate_faces=64, max_grid_voxels=8_000_000):
    """Sample an oriented, closest-triangle signed-distance field.

    Candidate faces are nearest by centroid, then evaluated with exact
    point-to-triangle distance. The finite candidate set is an approximation;
    spot-checks below reject meshes on which it misses closer triangles.
    """
    points = np.asarray(points, dtype=np.float64)
    triangles = points[faces]
    centers = triangles.mean(axis=1)
    directions = np.cross(triangles[:, 1] - triangles[:, 0],
                          triangles[:, 2] - triangles[:, 0])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    radius = np.max(np.linalg.norm(triangles - centers[:, None, :], axis=2), axis=1)
    maximum_radius = float(radius.max())
    margin = band_mm + 3 * spacing_mm
    origin = np.floor((points.min(axis=0) - margin) / spacing_mm) * spacing_mm
    upper = np.ceil((points.max(axis=0) + margin) / spacing_mm) * spacing_mm
    shape = np.rint((upper - origin) / spacing_mm).astype(int) + 1
    total = int(np.prod(shape))
    if total > max_grid_voxels:
        raise ValueError(f'Signed-distance grid needs {total:,} voxels; increase spacing or max_grid_voxels')
    if candidate_faces >= len(faces):
        candidate_faces = len(faces)
    centroid_tree = cKDTree(centers)
    vertex_tree = cKDTree(points)
    normals = vertex_normals(points, faces)
    output = np.empty(total, dtype=np.float32)
    near_count = 0
    sample_error = 0.
    sample_count = 0
    for start in range(0, total, 32768):
        stop = min(total, start + 32768)
        linear = np.arange(start, stop, dtype=np.int64)
        voxel = np.column_stack(np.unravel_index(linear, tuple(shape)))
        positions = origin + spacing_mm * voxel
        vertex_distance, vertex_id = vertex_tree.query(positions, workers=1)
        sign = np.einsum('ij,ij->i', positions - points[vertex_id], normals[vertex_id])
        block = np.where(sign >= 0, band_mm, -band_mm).astype(np.float32)
        nearby = vertex_distance <= band_mm + 2 * maximum_radius
        if nearby.any():
            nearest = np.asarray(centroid_tree.query(positions[nearby], k=candidate_faces,
                                                     workers=1)[1]).reshape(-1, candidate_faces)
            distances = _signed_distances(np.ascontiguousarray(positions[nearby]),
                                          np.ascontiguousarray(nearest), triangles, directions)
            block[nearby] = np.clip(distances, -band_mm, band_mm)
            near_count += int(nearby.sum())
            # Verify the finite centroid search against a wider candidate set.
            if candidate_faces < len(faces) and sample_count < 512:
                picked = np.linspace(0, len(distances)-1, min(8, len(distances)), dtype=int)
                wide_k = min(len(faces), max(4*candidate_faces, 256))
                wide_ids = np.asarray(centroid_tree.query(positions[nearby][picked],
                                                          k=wide_k, workers=1)[1]).reshape(-1, wide_k)
                wider = _signed_distances(np.ascontiguousarray(positions[nearby][picked]),
                                           np.ascontiguousarray(wide_ids), triangles, directions)
                sample_error = max(sample_error, float(np.max(np.abs(np.abs(distances[picked]) - np.abs(wider)))))
                sample_count += len(picked)
        output[start:stop] = block
    if sample_error > .05 * spacing_mm:
        raise ValueError(f'Closest-triangle candidate search missed by {sample_error:.3f} mm; increase candidate_faces')
    return output.reshape(tuple(shape)), origin, {'grid_shape': shape.tolist(),
        'grid_spacing_mm': float(spacing_mm), 'band_mm': float(band_mm),
        'candidate_faces': int(candidate_faces), 'near_voxels': near_count,
        'candidate_spotcheck_max_error_mm': sample_error}


def implicit_grid_principal_curvatures(points, faces, spacing_mm=1.,
                                        derivative_sigma_mm=1., band_mm=4.,
                                        candidate_faces=64, max_grid_voxels=8_000_000,
                                        gradient_fallback_threshold=.2):
    """Differentiate a mesh-derived signed distance grid at indexed vertices.

    Curvatures are the eigenvalues of one symmetric tangent-plane shape
    operator, so their mean and product are internally consistent.
    """
    if min(spacing_mm, derivative_sigma_mm, band_mm, gradient_fallback_threshold) <= 0:
        raise ValueError('Grid spacing, derivative sigma and band must be positive')
    if gradient_fallback_threshold > .5:
        raise ValueError('Gradient fallback threshold must be at most 0.5')
    if band_mm < 3 * derivative_sigma_mm:
        raise ValueError('Grid band must cover at least three derivative sigmas')
    grid, origin, details = _surface_signed_distance_grid(
        points, faces, float(spacing_mm), float(band_mm), int(candidate_faces),
        int(max_grid_voxels))
    query = ((np.asarray(points) - origin) / spacing_mm).T
    sigma = derivative_sigma_mm / spacing_mm
    gradient = np.empty((len(points), 3), dtype=np.float64)
    hessian = np.empty((len(points), 3, 3), dtype=np.float64)
    for a in range(3):
        orders = [0, 0, 0]
        orders[a] = 1
        derivative = ndimage.gaussian_filter(grid, sigma, order=orders,
                                             mode='nearest') / spacing_mm
        gradient[:, a] = ndimage.map_coordinates(derivative, query, order=1,
                                                  mode='nearest')
    for a in range(3):
        for b in range(a, 3):
            orders = [0, 0, 0]
            orders[a] += 1
            orders[b] += 1
            derivative = ndimage.gaussian_filter(grid, sigma, order=orders,
                                                 mode='nearest') / (spacing_mm**2)
            hessian[:, a, b] = ndimage.map_coordinates(derivative, query, order=1,
                                                        mode='nearest')
            hessian[:, b, a] = hessian[:, a, b]
    magnitude = np.linalg.norm(gradient, axis=1)
    if not np.isfinite(magnitude).all():
        raise ValueError('Signed-distance gradient is nonfinite at surface')
    # The grid remains directionally informative below magnitude 0.5 on the
    # cortical folds. Replacing all of those values with a local quadric fit
    # loses the grid's cross-bank information; reserve the fallback for nearly
    # flat gradients. Curvature amplification still has the original 0.5 floor.
    fallback = magnitude < gradient_fallback_threshold
    if np.mean(fallback) > .05:
        raise ValueError(f'Signed-distance gradient is weak at {np.mean(fallback):.1%} '
                         'of vertices; reduce derivative sigma or grid spacing')
    safe_magnitude = np.maximum(magnitude, .5)
    normals = gradient / np.maximum(magnitude[:, None], 1e-12)
    if fallback.any():
        normals[fallback] = vertex_normals(points, faces)[fallback]
    axis = np.zeros_like(normals)
    axis[:, 0] = 1.
    alternate = np.abs(normals[:, 0]) > .9
    axis[alternate] = (0., 1., 0.)
    tangent1 = np.cross(axis, normals)
    tangent1 /= np.linalg.norm(tangent1, axis=1)[:, None]
    tangent2 = np.cross(normals, tangent1)
    a = -np.einsum('ni,nij,nj->n', tangent1, hessian, tangent1) / safe_magnitude
    b = -np.einsum('ni,nij,nj->n', tangent1, hessian, tangent2) / safe_magnitude
    c = -np.einsum('ni,nij,nj->n', tangent2, hessian, tangent2) / safe_magnitude
    center = .5 * (a + c)
    difference = np.sqrt((.5 * (a-c))**2 + b*b)
    low, high = center-difference, center+difference
    if fallback.any():
        replacement_low, replacement_high, _ = quadric_principal_curvatures(
            points, faces, radius_mm=3., source_ids=np.flatnonzero(fallback))
        low[fallback] = replacement_low
        high[fallback] = replacement_high
    details.update(derivative_sigma_mm=float(derivative_sigma_mm),
                   gradient_magnitude_min=float(magnitude.min()),
                   gradient_magnitude_median=float(np.median(magnitude)),
                   weak_gradient_threshold=float(gradient_fallback_threshold),
                   quadric_fallback_vertices=int(fallback.sum()),
                   quadric_fallback_radius_mm=3.)
    if not np.isfinite(low).all() or not np.isfinite(high).all():
        raise ValueError('Nonfinite signed-distance curvature')
    return low, high, details
