import numpy as np
from scipy.spatial.transform import Rotation
from numpy.polynomial.legendre import leggauss
from swmtrack.harmonics import real_sh, rotate_to_local, polar_operator, circular_basis


def test_rotation_matches_direct_world_evaluation():
    rng = np.random.default_rng(8)
    for degree in (0, 4, 16):
        size = (degree+1)*(degree+2)//2
        coefficients = rng.normal(size=(7, size))
        frames = np.concatenate((Rotation.random(4, random_state=rng).as_matrix(),
            np.eye(3)[None], Rotation.from_euler("y", 180, degrees=True).as_matrix()[None],
            Rotation.from_euler("z", 60, degrees=True).as_matrix()[None]))
        points = rng.normal(size=(60, 3)); points /= np.linalg.norm(points, axis=1)[:, None]
        local = rotate_to_local(coefficients, frames, degree)
        for i in range(len(frames)):
            np.testing.assert_allclose(real_sh(points, degree)@local[i],
                real_sh(points@frames[i].T, degree)@coefficients[i], atol=2e-11)


def test_operator_matches_independent_unclipped_spherical_integral():
    rng = np.random.default_rng(19)
    degree = 16
    coefficients = rng.normal(size=153)
    angles = np.linspace(0, 2*np.pi, 47, endpoint=False)
    z, weights = leggauss(80)
    points = np.stack(np.broadcast_arrays(np.sqrt(1-z*z)[None, :]*np.cos(angles)[:, None],
        np.sqrt(1-z*z)[None, :]*np.sin(angles)[:, None], z[None, :]), axis=-1)
    explicit = (real_sh(points.reshape(-1, 3), degree)@coefficients).reshape(len(angles), -1)@weights
    direct = circular_basis(angles, degree)@polar_operator(degree)@coefficients
    np.testing.assert_allclose(direct, explicit, atol=2e-11)
    assert np.count_nonzero(polar_operator(16)) == 73


def test_constant_fod_and_antipodal_symmetry():
    c = np.zeros(153); c[0] = np.sqrt(4*np.pi)
    result = polar_operator(16)@c
    np.testing.assert_allclose(result, np.r_[2., np.zeros(16)], atol=1e-13)
    angles = np.linspace(-2, 3, 27)
    np.testing.assert_allclose(circular_basis(angles, 16), circular_basis(angles+np.pi, 16), atol=2e-14)
