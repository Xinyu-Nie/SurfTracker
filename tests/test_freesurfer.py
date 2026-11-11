import numpy as np
import nibabel as nib
import pytest
from nibabel.freesurfer.io import write_geometry, write_annot, write_morph_data
from swmtrack.freesurfer import import_freesurfer, surface_ras_to_scanner
from swmtrack.io import read_surface


def make_subject(tmp_path):
    root = tmp_path/"subject"
    for sub in ("mri", "surf", "label"):
        (root/sub).mkdir(parents=True)
    affine = np.diag([-1., 2., 3., 1.]); affine[:3, 3] = [-10, -20, -30]
    image = nib.MGHImage(np.zeros((20, 20, 20), dtype=np.float32), affine)
    nib.save(image, root/"mri/orig.mgz")
    scanner = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.]])
    source = nib.affines.apply_affine(np.linalg.inv(surface_ras_to_scanner(image)), scanner)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    volume_info = {"head": np.array([2, 0, 20]), "valid": "1", "filename": "orig.mgz",
                   "volume": np.array(image.shape), "voxelsize": nib.affines.voxel_sizes(image.affine),
                   "xras": [-1, 0, 0], "yras": [0, 1, 0], "zras": [0, 0, 1],
                   "cras": surface_ras_to_scanner(image)[:3, 3]}
    write_geometry(str(root/"surf/lh.white"), source, faces, volume_info=volume_info)
    write_geometry(str(root/"surf/lh.pial"), source+1, faces, volume_info=volume_info)
    table = np.array([[0, 0, 0, 0], [1, 2, 3, 0], [4, 5, 6, 0]])
    write_annot(str(root/"label/lh.aparc.annot"), np.array([1, 2, 1, 2]), table,
                [b"unknown", b"precentral", b"postcentral"])
    write_morph_data(str(root/"surf/lh.curv"), np.arange(4, dtype=np.float32))
    fod = tmp_path/"fod.nii"
    nib.save(nib.Nifti1Image(np.zeros((20, 20, 20, 1), dtype=np.float32), np.diag([1., 1., 1., 1.])), fod)
    return root, fod, scanner, faces


def test_import_preserves_physical_coordinates_and_atlas(tmp_path):
    root, fod, scanner, faces = make_subject(tmp_path)
    result = import_freesurfer(root, fod, tmp_path/"out", assume_scanner_aligned=True,
                               hemispheres=("lh",), gyrus_a="precentral", gyrus_b="postcentral")
    points, actual_faces, data = read_surface(tmp_path/"out/lh.white.world.vtk")
    np.testing.assert_allclose(points, scanner, atol=2e-6)
    np.testing.assert_array_equal(actual_faces, faces)
    np.testing.assert_array_equal(data["source_vertex_id"], np.arange(4))
    np.testing.assert_array_equal(data["gyrus_a"], [1, 0, 1, 0])
    np.testing.assert_array_equal(data["aparc"], [1+2*256+3*65536, 4+5*256+6*65536]*2)
    assert result["hemispheres"]["lh"]["surfaces"]["white"]["in_fod_domain_fraction"] == 1


def test_explicit_affine_and_reflection_recompute_outward_normals(tmp_path):
    root, fod, scanner, faces = make_subject(tmp_path)
    transform = np.diag([-1., 2., 1., 1.]); transform[:3, 3] = [10, 3, 5]
    import_freesurfer(root, fod, tmp_path/"out", scanner_to_world=transform, hemispheres=("lh",))
    points, actual_faces, data = read_surface(tmp_path/"out/lh.white.world.vtk")
    np.testing.assert_allclose(points, nib.affines.apply_affine(transform, scanner), atol=2e-6)
    np.testing.assert_array_equal(actual_faces, faces[:, [0, 2, 1]])
    assert np.all(np.einsum('ij,ij->i', points-points.mean(0), data["normals"]) > 0)


@pytest.mark.parametrize("mode", ["missing", "both", "singular"])
def test_import_requires_unambiguous_valid_transform(tmp_path, mode):
    root, fod, _, _ = make_subject(tmp_path)
    arguments = {} if mode == "missing" else (
        {"scanner_to_world": np.eye(4), "assume_scanner_aligned": True} if mode == "both" else
        {"scanner_to_world": np.zeros((4, 4))})
    with pytest.raises(ValueError, match="affine|Supply"):
        import_freesurfer(root, fod, tmp_path/"out", hemispheres=("lh",), **arguments)


def test_import_rejects_changed_pial_vertex_identity(tmp_path):
    root, fod, _, faces = make_subject(tmp_path)
    write_geometry(str(root/"surf/lh.pial"), np.ones((4, 3)), faces[::-1])
    with pytest.raises(ValueError, match="topology"):
        import_freesurfer(root, fod, tmp_path/"out", assume_scanner_aligned=True, hemispheres=("lh",))


def test_cli_import_uses_explicit_world_frame(tmp_path, capsys):
    import json
    from swmtrack.cli import main
    root, fod, scanner, _ = make_subject(tmp_path)
    main(["import-freesurfer", "--subject-dir", str(root), "--fod", str(fod),
          "--output-dir", str(tmp_path/"out"), "--hemisphere", "lh",
          "--assume-scanner-aligned", "--gyrus-a", "precentral"])
    result = json.loads(capsys.readouterr().out)
    assert result["hemispheres"]["lh"]["vertices"] == len(scanner)
