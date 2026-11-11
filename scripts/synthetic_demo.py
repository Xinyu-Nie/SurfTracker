"""Generate analytic inputs and run prepare -> track with no human subject data.

The flat surface and constant isotropic SH field demonstrate software I/O and
endpoint constraints. They do not model cortical anatomy or test segmentation.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from swmtrack.io import read_surface, write_surface, write_tractograms
from swmtrack.tracking import TrackingConfig, prepare, save_prepared, track
from swmtrack.validation import validate_exports


def check_tract_vtk(path, paths):
    """Read binary VTK points/connectivity independently of the writer."""
    with path.open('rb') as stream:
        assert stream.readline().startswith(b'# vtk DataFile Version')
        stream.readline()
        assert stream.readline().strip() == b'BINARY'
        assert stream.readline().strip() == b'DATASET POLYDATA'
        kind, n, dtype = stream.readline().split()
        assert kind == b'POINTS' and dtype == b'float'
        points = np.frombuffer(stream.read(int(n) * 12), dtype='>f4').reshape(-1, 3)
        np.testing.assert_array_equal(points, np.concatenate(paths).astype(np.float32))
        assert stream.readline() == b'\n'
        kind, count, size = stream.readline().split()
        assert kind == b'LINES' and int(count) == len(paths)
        cells = np.frombuffer(stream.read(int(size) * 4), dtype='>i4')
        offset = first = 0
        for points in paths:
            assert cells[first] == len(points)
            np.testing.assert_array_equal(cells[first+1:first+1+len(points)],
                                          np.arange(offset, offset+len(points)))
            first += len(points) + 1
            offset += len(points)
        assert first == len(cells) and stream.readline() == b'\n' and not stream.read()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--n-seeds', type=int, default=256)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be new or empty')
    config = TrackingConfig(n_seeds=args.n_seeds, threads=1)
    config.validate()
    output.mkdir(parents=True, exist_ok=True)

    n = 24
    x, y = np.meshgrid(np.linspace(2, 26, n + 1), np.linspace(2, 26, n + 1))
    points = np.column_stack((x.ravel(), y.ravel(), np.full(x.size, 4.0)))
    faces = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b, c = a + 1, a + n + 1
            faces.extend(((a, b, c + 1), (a, c + 1, c)))
    faces = np.asarray(faces, dtype=np.int32)
    left, right = points[:, 0] <= 5, points[:, 0] >= 23
    data = dict(roi=np.ones(len(points), np.int32),
                seed=((points[:, 0] >= 11) & (points[:, 0] <= 17)).astype(np.int32),
                gyrus_a=left.astype(np.int32), gyrus_b=right.astype(np.int32),
                gyrus=(left | right).astype(np.int32),
                sulcal_segmentation=(~(left | right)).astype(np.int32),
                source_vertex_id=np.arange(len(points), dtype=np.int32))
    surface, fod = output / 'surface.vtk', output / 'fod.nii.gz'
    write_surface(surface, points, faces, data)
    # l=0 coefficient sqrt(4*pi) gives the constant spherical field F(n)=1.
    image = nib.Nifti1Image(np.full((30, 30, 10, 1), np.sqrt(4*np.pi), np.float32), np.eye(4))
    nib.save(image, fod)
    actual_points, actual_faces, actual_data = read_surface(surface)
    np.testing.assert_array_equal(actual_points, points)
    np.testing.assert_array_equal(actual_faces, faces)
    for key in data:
        np.testing.assert_array_equal(actual_data[key], data[key])

    prepared = prepare(surface, fod, roi_output=output / 'tracking_surface.vtk')
    save_prepared(output / 'prepared.npz', prepared)
    paths, report, diagnostics = track(prepared, config)
    if not paths:
        raise RuntimeError('Synthetic example unexpectedly produced no accepted tracks')
    accepted = diagnostics['accepted_trial_ids']
    np.testing.assert_array_equal(np.sort(diagnostics['endpoint_labels'][accepted], axis=1),
                                  np.tile([1, 2], (len(accepted), 1)))
    if not np.all((diagnostics['lengths_mm'][accepted] >= config.min_length_mm) &
                  (diagnostics['lengths_mm'][accepted] <= config.max_length_mm)):
        raise RuntimeError('Accepted tracks violate configured length limits')
    outputs = write_tractograms(output / 'tracks', paths, image, vtk=True)
    validation = validate_exports(output / 'tracks', paths, image, prepared, diagnostics)
    check_tract_vtk(output / 'tracks.vtk', paths)
    validation['vtk_points_and_polylines_exact'] = True
    result = dict(status='passed', human_subject_data=False, sh_basis='legacy-real-even',
                  scope='analytic prepare/track/export smoke example; no FreeSurfer import or segmentation',
                  tracking_config=asdict(config), tracking=report, validation=validation, outputs=outputs)
    (output / 'demo.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='passed', accepted=len(paths), output_dir=str(output)), indent=2))


if __name__ == '__main__':
    main()
