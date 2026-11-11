"""Portable command-line helpers require only this source archive and runtime."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import numpy as np
import pytest
from nibabel.freesurfer.io import write_geometry
import swmtrack
from swmtrack.io import read_surface, write_surface, sha256
from test_preprocessing import octahedron


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def run_helper(name, *args):
    env = dict(os.environ, PYTHONPATH=str(Path(swmtrack.__file__).resolve().parents[1]),
               OPENBLAS_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="swmtrack-helper-") as workdir:
        return subprocess.run([sys.executable, str(SCRIPTS / name), *map(str, args)],
                              cwd=workdir, env=env, capture_output=True, text=True)


@pytest.mark.parametrize("merged", [False, True])
def test_qc_exports_original_labels_and_inflation_without_modifying_run(tmp_path, merged):
    run = tmp_path / "run"
    pre = run / "lh/preprocessing"
    pre.mkdir(parents=True)
    (run / "surfaces").mkdir()
    subject = tmp_path / "subject"
    (subject / "surf").mkdir(parents=True)
    points, faces = octahedron()
    atlas = np.array([0, 17, 17, 42, 42, 0], np.int32)
    segmentation = np.array([0, 0, 1, 1, 0, 0], np.int32)
    roi = atlas > 0
    a = (segmentation == 0) & (roi if merged else atlas == 17)
    b = np.zeros(6, bool) if merged else (segmentation == 0) & (atlas == 42)
    data = {"aparc": atlas, "roi": roi, "seed": roi & (segmentation == 1),
            "sulcal_segmentation": segmentation, "gyrus_a": a, "gyrus_b": b,
            "source_vertex_id": np.arange(6, dtype=np.int32)}
    source = pre / "preprocessed.world.vtk"
    write_surface(source, points, faces, data)
    original = sha256(source)
    write_geometry(str(subject / "surf/lh.inflated"), (points * 2).astype(np.float32), faces)
    (run / "run.json").write_text(json.dumps({"status": "complete", "hemispheres": {"lh": {}},
        "signature": {"annotation": "aparc", "subject_dir": str(subject)}}))
    (run / "surfaces/freesurfer_import.json").write_text(json.dumps({"hemispheres": {"lh": {
        "labels": [{"id": 17, "rgba": [20, 30, 40, 255]}, {"id": 42, "rgba": [50, 60, 70, 255]}]}}}))
    out = tmp_path / "qc"
    result = run_helper("export_qc.py", "--run-dir", run, "--output-dir", out, "--inflated")
    assert result.returncode == 0, result.stderr
    assert sha256(source) == original
    for kind, expected_points in [("white", points), ("inflated", points * 2)]:
        p, f, d = read_surface(out / "lh" / f"labels.{kind}.vtk")
        np.testing.assert_array_equal(p, expected_points)
        np.testing.assert_array_equal(f, faces)
        for key, value in data.items():
            np.testing.assert_array_equal(d[key], value)
        np.testing.assert_array_equal(d["candidate_segmentation"], [2, 0, 1, 1, 0, 2])
        np.testing.assert_array_equal(d["endpoint_region"], a.astype(int) + 2 * b)
        np.testing.assert_array_equal(d["anatomical_roi"], atlas)
    repeat = run_helper("export_qc.py", "--run-dir", run, "--output-dir", out)
    assert repeat.returncode != 0 and "nonempty" in repeat.stderr


@pytest.mark.parametrize("bad_path", [None, "swmtrack/_native/bin/GeoDist", "swmtrack/private.nii.gz", "cohort/helper.py"])
def test_archive_guard_accepts_source_only_and_rejects_unwanted_payload(tmp_path, bad_path):
    wheel = tmp_path / "sample.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("swmtrack/__init__.py", '__version__ = "0.3.0.dev12"\n')
        z.writestr("surftracker-0.3.0.dev12.dist-info/METADATA", "Metadata-Version: 2.4\nName: SurfTracker\nVersion: 0.3.0.dev12\n\n")
        if bad_path:
            z.writestr(bad_path, b"not a permitted release file")
    result = run_helper("check_distributions.py", wheel)
    if bad_path:
        assert result.returncode != 0
    else:
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)["archives"][0]["project_license_ready"] is False
        result = run_helper("check_distributions.py", wheel, "--require-license")
        assert result.returncode != 0 and "license" in result.stderr
