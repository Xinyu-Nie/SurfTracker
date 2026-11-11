import numpy as np
import pytest
from nibabel.affines import apply_affine
from swmtrack.registration import read_itk_affine_lps, itk_resampling_affine_to_ras_points


def write_transform(path, matrix, translation, center):
    path.write_text("#Insight Transform File V1.0\n#Transform 0\n"
                    "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "+
                    " ".join(map(str, [*matrix.ravel(), *translation]))+
                    "\nFixedParameters: "+" ".join(map(str, center))+"\n")


def test_centered_itk_affine_and_image_point_direction(tmp_path):
    path = tmp_path/"centered.txt"
    matrix = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    center, translation = np.array([10., 20., 30.]), np.array([2., -3., 4.])
    write_transform(path, matrix, translation, center)
    x = np.array([[5., 7., 9.], [11., 20., 35.]])
    moving = (x-center)@matrix.T+center+translation
    np.testing.assert_allclose(apply_affine(read_itk_affine_lps(path), x), moving)
    flip = np.array([-1., -1., 1.])
    np.testing.assert_allclose(apply_affine(itk_resampling_affine_to_ras_points(path), moving*flip), x*flip)


def test_itk_translation_requires_inverse_and_lps_ras_flip(tmp_path):
    path = tmp_path/"translation.txt"
    write_transform(path, np.eye(3), [2., -3., 1.], [0., 0., 0.])
    np.testing.assert_allclose(apply_affine(itk_resampling_affine_to_ras_points(path), [[5., 6., 7.]]), [[7., 3., 6.]])


@pytest.mark.parametrize("text", ["", "Transform: BSplineTransform_double_3_3\n",
                                    "Transform: AffineTransform_double_3_3\nParameters: 1 2\n",
                                    "Transform: AffineTransform_double_3_3\nTransform: AffineTransform_double_3_3\n"])
def test_reject_nonaffine_or_incomplete_itk_transform(tmp_path, text):
    path = tmp_path/"bad.txt";path.write_text(text)
    with pytest.raises(ValueError):
        read_itk_affine_lps(path)
