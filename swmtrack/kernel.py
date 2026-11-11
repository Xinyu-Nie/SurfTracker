"""Compiled, deterministic propagation with fixed millimeter integration steps."""

import numpy as np
from numba import njit, prange

STOP_NAMES = {0: "not_started", 1: "gyrus_a", 2: "gyrus_b", 3: "boundary", 4: "no_supported_direction",
              5: "max_length", 6: "point_guard", 7: "vertex_hit", 8: "invalid_geometry", 9: "overlapping_targets"}


@njit(cache=True, inline="always")
def local_xy(difference, frame):
    x = difference[0]*frame[0, 0]+difference[1]*frame[1, 0]+difference[2]*frame[2, 0]
    y = difference[0]*frame[0, 1]+difference[1]*frame[1, 1]+difference[2]*frame[2, 1]
    return x, y


@njit(cache=True)
def random_uniform(state):
    """SplitMix64: explicit per-trial state, independent of worker scheduling."""
    state = state + np.uint64(0x9E3779B97F4A7C15)
    z = state
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z ^= z >> np.uint64(31)
    return state, float(z >> np.uint64(11)) * (1.0 / 9007199254740992.0)


@njit(cache=True)
def categorical(weights, state):
    total = 0.0
    for value in weights:
        if np.isfinite(value) and value > 0:
            total += value
    if total <= 0 or not np.isfinite(total):
        return -1, state
    state, u = random_uniform(state)
    limit, cumulative, last = u*total, 0.0, -1
    for i in range(len(weights)):
        if np.isfinite(weights[i]) and weights[i] > 0:
            cumulative += weights[i]
            last = i
            if limit < cumulative:
                return i, state
    return last, state


@njit(cache=True)
def ray_exit(xy, x, y, dx, dy):
    """Direct 2D ray/edge intersection. A vertex tie is deterministic and flagged."""
    best, edge, vertex = np.inf, -1, False
    for e in range(3):
        j = (e+1) % 3
        ax, ay = xy[e, 0], xy[e, 1]
        rx, ry = xy[j, 0]-ax, xy[j, 1]-ay
        det = dx*ry-dy*rx
        if abs(det) < 1e-12*np.sqrt(rx*rx+ry*ry):
            continue
        qx, qy = ax-x, ay-y
        t = (qx*ry-qy*rx)/det
        u = (qx*dy-qy*dx)/det
        if t > 1e-9 and u >= -1e-10 and u <= 1+1e-10 and t < best-1e-10:
            best, edge = t, e
            vertex = u < 1e-10 or u > 1-1e-10
    return best, edge, vertex


@njit(cache=True)
def interpolate_coefficients(coefficients, xy, x, y):
    out = np.empty(coefficients.shape[1])
    if coefficients.shape[0] == 1:
        for k in range(len(out)):
            out[k] = coefficients[0, k]
    else:
        b2 = y/xy[2, 1]
        b1 = (x-xy[2, 0]*b2)/xy[1, 0]
        b0 = 1-b1-b2
        for k in range(len(out)):
            out[k] = b0*coefficients[0, k]+b1*coefficients[1, k]+b2*coefficients[2, k]
    return out


@njit(cache=True)
def circular_value(coefficients, angle):
    value = coefficients[0]
    cosine, sine = np.cos(2*angle), np.sin(2*angle)
    cm, sm = cosine, sine
    for k in range((len(coefficients)-1)//2):
        value += coefficients[2*k+1]*cm+coefficients[2*k+2]*sm
        cm, sm = cm*cosine-sm*sine, sm*cosine+cm*sine
    return max(value, 0.0)  # Deliberate clipping AFTER the 2D projection.


@njit(cache=True)
def sample_direction(coefficients, xy, x, y, center, half_width, resolution, min_fod, state):
    count = max(1, int(np.ceil(2*half_width/resolution)))
    width = 2*half_width/count
    weights = np.zeros(count)
    local = interpolate_coefficients(coefficients, xy, x, y)
    for k in range(count):
        angle = center-half_width+(k+0.5)*width
        value = circular_value(local, angle)
        if not np.isfinite(value) or value < min_fod or value <= 0:
            continue
        distance, edge, _ = ray_exit(xy, x, y, np.cos(angle), np.sin(angle))
        if edge >= 0:
            weights[k] = value
    chosen, state = categorical(weights, state)
    if chosen < 0:
        return 0.0, 0.0, False, state
    angle = center-half_width+(chosen+0.5)*width
    return np.cos(angle), np.sin(angle), True, state


@njit(cache=True)
def trace_half(points, faces, frames, xy, neighbors, transport, coefficients, targets,
               seed_face, seed_point, dx, dy, step_mm, turn_radians, angular_resolution,
               min_fod, max_length, max_points, stop_at_targets, state):
    path = np.empty((max_points, 3))
    path[0] = seed_point
    count, face, length, samples, crossings = 1, seed_face, 0.0, 0, 0
    origin = points[faces[face, 0]]
    difference = seed_point-origin
    x, y = local_xy(difference, frames[face])
    remaining = min(step_mm, max_length)
    status, end_face = 6, face
    while count < max_points:
        if max_length-length <= 1e-9:
            status = 5
            break
        if remaining <= 1e-9:
            dx, dy, valid, state = sample_direction(coefficients[face], xy[face], x, y,
                np.arctan2(dy, dx), turn_radians, angular_resolution, min_fod, state)
            samples += 1
            if not valid:
                status = 4
                break
            remaining = min(step_mm, max_length-length)
        distance, edge, vertex = ray_exit(xy[face], x, y, dx, dy)
        if edge < 0 or not np.isfinite(distance):
            status = 8
            break
        travel = min(distance, remaining, max_length-length)
        x, y = x+travel*dx, y+travel*dy
        origin = points[faces[face, 0]]
        for dim in range(3):
            path[count, dim] = origin[dim]+x*frames[face, dim, 0]+y*frames[face, dim, 1]
        count += 1
        length += travel
        remaining -= travel
        end_face = face
        if distance <= travel+1e-9:
            crossings += 1
            next_face = neighbors[face, edge]
            if vertex:
                status = 7
                break
            if next_face < 0:
                status = 3
                break
            if stop_at_targets and targets[next_face] > 0:
                label = targets[next_face]
                status = label if label < 3 else 9
                end_face = next_face
                break
            tx = transport[face, edge, 0, 0]*dx+transport[face, edge, 0, 1]*dy
            ty = transport[face, edge, 1, 0]*dx+transport[face, edge, 1, 1]*dy
            norm = np.sqrt(tx*tx+ty*ty)
            if norm < 1e-12:
                status = 8
                break
            dx, dy = tx/norm, ty/norm
            face = next_face
            difference = path[count-1]-points[faces[face, 0]]
            x, y = local_xy(difference, frames[face])
    return path, count, status, end_face, length, state, samples, crossings


@njit(cache=True, parallel=True)
def track_batch(points, faces, frames, xy, neighbors, transport, coefficients, targets,
                seed_faces, seed_points, trial_ids, global_seed, step_mm, turn_radians,
                angular_resolution, min_fod, max_length, max_points, stop_at_targets):
    batch = len(seed_faces)
    paths = np.empty((batch, 2*max_points-1, 3))
    counts = np.zeros(batch, dtype=np.int64)
    stops = np.zeros((batch, 2), dtype=np.int64)
    end_faces = np.full((batch, 2), -1, dtype=np.int64)
    lengths = np.zeros(batch)
    operations = np.zeros((batch, 2), dtype=np.int64)
    for i in prange(batch):
        state = np.uint64(global_seed) ^ ((np.uint64(trial_ids[i])+np.uint64(1))*np.uint64(0xD2B74407B1CE6E93))
        face = seed_faces[i]
        difference = seed_points[i]-points[faces[face, 0]]
        x, y = local_xy(difference, frames[face])
        dx, dy, valid, state = sample_direction(coefficients[face], xy[face], x, y, 0.0,
            np.pi, angular_resolution, min_fod, state)
        if not valid:
            stops[i, 0], stops[i, 1] = 4, 4
            continue
        p, np_, sp, ep, lp, state, samples_p, crossings_p = trace_half(points, faces, frames, xy,
            neighbors, transport, coefficients, targets, face, seed_points[i], dx, dy, step_mm,
            turn_radians, angular_resolution, min_fod, max_length, max_points, stop_at_targets, state)
        n, nn, sn, en, ln, state, samples_n, crossings_n = trace_half(points, faces, frames, xy,
            neighbors, transport, coefficients, targets, face, seed_points[i], -dx, -dy, step_mm,
            turn_radians, angular_resolution, min_fod, max_length, max_points, stop_at_targets, state)
        for j in range(nn):
            paths[i, j] = n[nn-1-j]
        for j in range(1, np_):
            paths[i, nn+j-1] = p[j]
        counts[i] = nn+np_-1
        stops[i, 0], stops[i, 1] = sn, sp
        end_faces[i, 0], end_faces[i, 1] = en, ep
        lengths[i] = ln+lp
        operations[i, 0], operations[i, 1] = samples_n+samples_p+1, crossings_n+crossings_p
    return paths, counts, stops, end_faces, lengths, operations
