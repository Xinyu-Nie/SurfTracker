"""Verify an installed binary package and its synthetic example report."""
import argparse
from dataclasses import fields
from importlib.metadata import distribution
import inspect
import json
from pathlib import Path

import numpy as np
import swmtrack
from swmtrack import kernel, tracking, preprocessing, geometry, harmonics, workflow
from swmtrack._provenance import build_info


def validate(demo, version):
    package = distribution('SurfTracker')
    assert package.version == swmtrack.__version__ == version
    entry = next(e for e in package.entry_points if e.group == 'console_scripts' and e.name == 'swmtrack')
    assert entry.value == 'swmtrack.cli:main'
    for module in (kernel, tracking, preprocessing, geometry, harmonics, workflow):
        assert module.__file__.endswith(('.so', '.pyd')), module.__file__
        assert not Path(module.__file__).with_name(module.__name__.split('.')[-1]+'.py').exists()
    assert build_info()['implementation'] == 'cython-native' and kernel.NATIVE
    assert next(f for f in fields(tracking.TrackingConfig) if f.name == 'seed').type is int
    assert next(f for f in fields(preprocessing.PreprocessConfig) if f.name == 'smoothing_mass').type is str
    assert next(f for f in fields(tracking.Prepared) if f.name == 'coefficients').type is np.ndarray
    assert 'whole_hemisphere' in inspect.signature(workflow.run_pipeline).parameters
    config = preprocessing.PreprocessConfig.recommended()
    assert (config.smoothing_mass, config.segmentation_energy) == ('consistent', 'robust-unit')
    assert config.between_gyri_seed_policy == 'connected-intersection'
    assert config.within_gyrus_seed_policy == 'intersection'
    result = json.loads(Path(demo).read_text())
    assert result['status'] == 'passed'
    assert result['tracking']['attempts'] == 256 and result['tracking']['accepted'] == 98
    assert result['validation']['streamline_digest'] == '8c4b8046e76d972584759a17a05cdef0f5d929214c2a538b783c4253103266bb'
    assert result['validation']['vtk_points_and_polylines_exact']
    return {'status': 'passed', 'version': version, 'command': entry.name,
            'core_implementation': 'native', 'synthetic_attempts': 256, 'synthetic_accepted': 98,
            'streamline_digest': result['validation']['streamline_digest']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--demo', type=Path, required=True)
    parser.add_argument('--version', default='0.3.0.dev18')
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    result = json.dumps(validate(args.demo, args.version), indent=2) + '\n'
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(result)
    print(result, end='')
