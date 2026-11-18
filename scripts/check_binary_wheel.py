"""Fail-closed payload and RECORD audit for compiled SurfTracker wheels."""
import argparse
import base64
import csv
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import zipfile

CORE = {'curvature_estimators', 'geometry', 'harmonics', 'kernel', 'open_preprocessing',
        'parcellation', 'preprocessing', 'tracking', 'workflow', '_segmentation_kernels'}
PUBLIC = {'__init__', '__main__', 'cli', 'io', 'freesurfer', 'registration',
          'validation', '_ants_registration', '_provenance', '_api'}
NOTICES = {'AUTHORS.rst', 'THIRD_PARTY_NOTICES.rst', 'CITATION.cff', 'CITATIONS.bib'}


def check(wheel):
    wheel = Path(wheel)
    if wheel.suffix != '.whl' or wheel.name.endswith('-none-any.whl'):
        raise ValueError('Expected a platform-specific compiled wheel')
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError('Duplicate wheel entries')
        files = {name: archive.read(name) for name in names if not name.endswith('/')}
    metadata_names = [name for name in files if name.endswith('.dist-info/METADATA')]
    if len(metadata_names) != 1:
        raise ValueError('Expected one distribution metadata directory')
    metadata_path = PurePosixPath(metadata_names[0]).parent
    metadata = BytesParser().parsebytes(files[metadata_names[0]])
    if metadata['Name'] != 'SurfTracker':
        raise ValueError('Unexpected distribution name')
    native = {}
    readable = set()
    for name, data in files.items():
        path = PurePosixPath(name)
        parts = path.parts
        if path.is_absolute() or '..' in parts or '\\' in name:
            raise ValueError('Unsafe wheel path: ' + name)
        if len(parts) >= 2 and parts[0] == str(metadata_path):
            if len(parts) == 2 and parts[1] in {'METADATA', 'WHEEL', 'RECORD', 'entry_points.txt', 'top_level.txt'}:
                continue
            if len(parts) == 3 and parts[1] == 'licenses' and parts[2] in NOTICES:
                continue
            raise ValueError('Unexpected metadata payload: ' + name)
        if len(parts) == 3 and parts[:2] == ('swmtrack', '_notices') and parts[2] in NOTICES:
            continue
        if len(parts) != 2 or parts[0] != 'swmtrack':
            raise ValueError('Unexpected wheel payload: ' + name)
        if parts[1] == '_build_manifest.json':
            continue
        if path.suffix == '.py' and path.stem in PUBLIC:
            readable.add(path.stem)
            continue
        module = parts[1].split('.')[0]
        if module in CORE and path.suffix in ('.so', '.pyd') and module not in native:
            if data[:4] != b'\x7fELF' and data[:2] != b'MZ' and data[:4] not in (b'\xcf\xfa\xed\xfe', b'\xca\xfe\xba\xbe'):
                raise ValueError('Invalid native module header: ' + name)
            # Exclude host paths and embedded source/code payloads, including
            # mistakes introduced by changed compiler/debugging options.
            if re.search(rb'/(?:ifs|ifshome|home|Users|tmp)/[ -~]{3,}', data):
                raise ValueError('Private absolute build path in binary: ' + name)
            if any(marker in data for marker in (b'def trace_half(', b'@njit(', b'def robust_unit_costs(')):
                raise ValueError('Embedded protected source in binary: ' + name)
            native[module] = hashlib.sha256(data).hexdigest()
            continue
        raise ValueError('Unapproved source, binary or data payload: ' + name)
    if set(native) != CORE or readable != PUBLIC:
        raise ValueError('Missing expected native modules or public interfaces')
    manifest = json.loads(files['swmtrack/_build_manifest.json'])
    if manifest.get('implementation') != 'cython-native' or set(manifest['compiled_modules']) != CORE:
        raise ValueError('Invalid compiled build manifest')
    if b'Root-Is-Purelib: false' not in files[str(metadata_path / 'WHEEL')]:
        raise ValueError('Native wheel declared pure Python')
    if b'swmtrack = swmtrack.cli:main' not in files[str(metadata_path / 'entry_points.txt')]:
        raise ValueError('swmtrack command is missing or changed')
    record_path = str(metadata_path / 'RECORD')
    rows = list(csv.reader(io.StringIO(files[record_path].decode())))
    if {row[0] for row in rows} != set(files) or len(rows) != len(files):
        raise ValueError('RECORD entries do not match the payload')
    for name, digest, size in rows:
        if name == record_path:
            if digest or size:
                raise ValueError('RECORD must not hash itself')
            continue
        expected = base64.urlsafe_b64encode(hashlib.sha256(files[name]).digest()).decode().rstrip('=')
        if digest != 'sha256=' + expected or int(size) != len(files[name]):
            raise ValueError('RECORD checksum mismatch: ' + name)
    return {'status': 'passed', 'wheel': str(wheel.resolve()), 'version': metadata['Version'],
            'sha256': hashlib.sha256(wheel.read_bytes()).hexdigest(), 'bytes': wheel.stat().st_size,
            'native_sha256': native, 'public_python_modules': sorted(readable),
            'protected_source_files': 0, 'files': len(files), 'build': manifest}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    output = json.dumps(check(args.wheel), indent=2) + '\n'
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
