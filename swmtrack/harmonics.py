"""Real SH rotations and the exact linear polar-marginal coefficient map.

SH convention: even l, m=-l..l; negative m is sqrt(2) Im(Y_l,|m|),
positive m is sqrt(2) Re(Y_l,m), including the Condon--Shortley phase.
This is the basis used by the supplied MATLAB matrixB.m. No clipping occurs
in 3D. Tracking evaluates max(Fourier(theta), 0) after projection to 2D.
"""

from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln, lpmv


def sh_indices(degree):
    if degree < 0 or degree % 2:
        raise ValueError("SH degree must be a nonnegative even integer")
    return [(ell, m) for ell in range(0, degree + 1, 2) for m in range(-ell, ell + 1)]


def degree_from_count(count):
    degree = int(round((np.sqrt(8 * count + 1) - 3) / 2))
    if degree % 2 or (degree + 1) * (degree + 2) // 2 != count:
        raise ValueError(f"{count} is not an even-degree SH coefficient count")
    return degree


def real_sh(directions, degree):
    directions = np.asarray(directions, dtype=float)
    directions = directions / np.linalg.norm(directions, axis=1)[:, None]
    z = np.clip(directions[:, 2], -1, 1)
    azimuth = np.arctan2(directions[:, 1], directions[:, 0])
    out = np.empty((len(directions), len(sh_indices(degree))))
    for j, (ell, m) in enumerate(sh_indices(degree)):
        k = abs(m)
        scale = np.sqrt((2 * ell + 1) / (4 * np.pi))
        if k:
            scale *= np.sqrt(2) * np.exp(0.5 * (gammaln(ell-k+1) - gammaln(ell+k+1)))
        angular = 1 if not k else (np.sin(k*azimuth) if m < 0 else np.cos(k*azimuth))
        out[:, j] = scale * lpmv(k, ell, z) * angular
    return out


def circular_basis(angles, degree):
    angles = np.atleast_1d(angles)
    out = np.ones((len(angles), degree + 1))
    for k in range(1, degree // 2 + 1):
        out[:, 2*k-1] = np.cos(2*k*angles)
        out[:, 2*k] = np.sin(2*k*angles)
    return out


@lru_cache(maxsize=8)
def polar_operator(degree):
    """Map local real SH to [DC, cos(2t), sin(2t), ...], without fitting.

    Even azimuthal orders have polynomial polar integrands of degree <= L;
    Gauss--Legendre integration is exact. Odd orders cancel by symmetry.
    """
    z, weights = leggauss(degree // 2 + 1)
    out = np.zeros((degree + 1, len(sh_indices(degree))))
    for j, (ell, m) in enumerate(sh_indices(degree)):
        if m % 2:
            continue
        k = abs(m)
        scale = np.sqrt((2*ell+1)/(4*np.pi))
        if k:
            scale *= np.sqrt(2) * np.exp(0.5*(gammaln(ell-k+1)-gammaln(ell+k+1)))
        row = 0 if m == 0 else k - 1 + int(m < 0)
        out[row, j] = scale * np.dot(weights, lpmv(k, ell, z))
    out[np.abs(out) < 1e-13] = 0
    out.setflags(write=False)
    return out


@lru_cache(maxsize=8)
def y90_blocks(degree):
    """Real-SH rotation blocks for f(R_y(pi/2) u), computed once.

    Exact spherical quadrature avoids per-face Euler/Wigner factorial sums.
    Products have polynomial degree <= 2L and azimuthal frequency <= 2L.
    """
    z, w = leggauss(degree + 1)
    theta = np.arange(2*degree+1) * (2*np.pi/(2*degree+1))
    zz, tt = np.meshgrid(z, theta, indexing="ij")
    radial = np.sqrt(1-zz*zz)
    points = np.column_stack((radial.ravel()*np.cos(tt.ravel()), radial.ravel()*np.sin(tt.ravel()), zz.ravel()))
    rotated = points[:, [2, 1, 0]].copy()
    rotated[:, 2] *= -1
    a, b = real_sh(points, degree), real_sh(rotated, degree)
    weights = np.repeat(w, len(theta)) * (2*np.pi/len(theta))
    blocks = []
    offset = 0
    for ell in range(0, degree+1, 2):
        size = 2*ell+1
        block = a[:, offset:offset+size].T @ (weights[:, None]*b[:, offset:offset+size])
        if not np.allclose(block.T @ block, np.eye(size), atol=2e-11):
            raise RuntimeError("SH rotation quadrature lost orthogonality")
        blocks.append(block)
        offset += size
    return tuple(blocks)


def _z_rotate(values, angle, ell):
    out = values.copy()
    for m in range(1, ell+1):
        cosine, sine = np.cos(m*angle), np.sin(m*angle)
        neg, pos = values[:, ell-m], values[:, ell+m]
        out[:, ell+m] = pos*cosine + neg*sine
        out[:, ell-m] = neg*cosine - pos*sine
    return out


def rotate_to_local(coefficients, frames, degree):
    """frames[n,:,j] is local axis j in the SH coefficient direction frame."""
    coefficients = np.asarray(coefficients, dtype=float)
    frames = np.asarray(frames, dtype=float)
    if coefficients.shape != (len(frames), len(sh_indices(degree))):
        raise ValueError("Coefficient/frame shape mismatch")
    beta = np.arccos(np.clip(frames[:, 2, 2], -1, 1))
    alpha = np.arctan2(frames[:, 1, 2], frames[:, 0, 2])
    gamma = np.arctan2(frames[:, 2, 1], -frames[:, 2, 0])
    north = beta < 1e-10
    south = np.pi-beta < 1e-10
    alpha[north] = np.arctan2(frames[north, 1, 0], frames[north, 0, 0])
    alpha[south] = np.arctan2(-frames[south, 1, 0], -frames[south, 0, 0])
    gamma[north | south] = 0
    out = np.empty_like(coefficients)
    offset = 0
    for ell, y90 in zip(range(0, degree+1, 2), y90_blocks(degree)):
        size = 2*ell+1
        block = coefficients[:, offset:offset+size]
        block = _z_rotate(block, alpha-np.pi/2, ell)
        block = block @ y90  # T_y(-pi/2), row-vector convention
        block = _z_rotate(block, beta, ell)
        block = block @ y90.T
        block = _z_rotate(block, gamma+np.pi/2, ell)
        out[:, offset:offset+size] = block
        offset += size
    return out


def project_coefficients(coefficients, frames, degree):
    return np.ascontiguousarray(rotate_to_local(coefficients, frames, degree) @ polar_operator(degree).T)
