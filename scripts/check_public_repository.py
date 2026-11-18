"""Check the public documentation/download tree without importing the core."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
from urllib.parse import unquote, urlsplit
import zipfile

from check_binary_wheel import check as check_wheel

PUBLIC_HELPERS = {'scripts/check_binary_wheel.py', 'scripts/check_public_repository.py',
                  'scripts/validate_installation.py', 'scripts/export_qc.py',
                  'scripts/synthetic_demo.py'}
ROOT_FILES = {'.gitignore', 'README.md', 'README.rst', 'AUTHORS.rst', 'CHANGELOG.rst',
              'CITATION.cff', 'CITATIONS.bib', 'THIRD_PARTY_NOTICES.rst'}
DOWNLOAD_FILES = {'manifest.json', 'SHA256SUMS', 'requirements-validation.txt'}
CORE_DIRECTORIES = {'swmtrack', 'build_support', 'tests', 'research', 'research_tests',
                    'Surface_tracking_c++', 'Surface_tracking_matlab'}


def check(root, check_rst=False):
    root = Path(root).resolve()
    # Inspect the candidate working tree as well as new, unignored files. Git
    # history is deliberately outside this current-tree audit.
    paths = subprocess.check_output(['git', '-C', str(root), 'ls-files',
        '--cached', '--others', '--exclude-standard', '-z']).decode().split('\0')
    paths = sorted({name for name in paths if name and (root/name).is_file()})
    for name in paths:
        p = PurePosixPath(name)
        parts = p.parts
        if p.is_absolute() or '..' in parts or parts[0] in CORE_DIRECTORIES:
            raise ValueError('Protected/unexpected path: ' + name)
        allowed = (name in ROOT_FILES or name in PUBLIC_HELPERS
            or name == '.github/workflows/ci.yml'
            or parts[0] in ('docs', 'examples') and p.suffix in ('.md', '.rst')
            or name == 'binaries/README.md'
            or len(parts) == 3 and parts[0] == 'binaries'
                and (p.suffix in ('.whl', '.zip') or parts[2] in DOWNLOAD_FILES))
        if not allowed:
            raise ValueError('Unapproved public payload: ' + name)
        if p.suffix in ('.md', '.rst', '.yml', '.cff'):
            text = (root/name).read_text()
            if re.search(r'/(?:ifs|ifshome)/', text):
                raise ValueError('Institutional path in public document: ' + name)
            if p.suffix == '.md':
                for target in re.findall(r'\]\(([^)]+)\)', text):
                    target = target.strip('<>')
                    parsed = urlsplit(target)
                    if parsed.scheme or parsed.netloc or not parsed.path:
                        continue
                    destination = (root/name).parent / unquote(parsed.path)
                    if not destination.exists():
                        raise ValueError(f'Broken local link in {name}: {target}')
    versions = []
    kit_rst_documents = 0
    for path in sorted((root/'binaries').glob('*/manifest.json')):
        manifest = json.loads(path.read_text())
        directory = path.parent
        names = set(manifest['files'])
        actual = {p.name for p in directory.iterdir() if p.is_file()}
        if actual != names | {'manifest.json', 'SHA256SUMS'}:
            raise ValueError('Download directory differs from manifest')
        for name, record in manifest['files'].items():
            data = (directory/name).read_bytes()
            if hashlib.sha256(data).hexdigest() != record['sha256'] or len(data) != record['bytes']:
                raise ValueError('Download hash/size mismatch: ' + name)
        for line in (directory/'SHA256SUMS').read_text().splitlines():
            digest, name = line.split('  ', 1)
            if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
                raise ValueError('SHA256SUMS mismatch: ' + name)
        wheels = list(directory.glob('*.whl'))
        if len(wheels) != 1:
            raise ValueError('Expected one wheel per supported binary version')
        wheel_report = check_wheel(wheels[0])
        if wheel_report['version'] != manifest['version']:
            raise ValueError('Wheel/manifest versions differ')
        kits = list(directory.glob('*.zip'))
        if len(kits) != 1:
            raise ValueError('Expected one complete binary kit')
        with zipfile.ZipFile(kits[0]) as archive:
            for name in archive.namelist():
                p = PurePosixPath(name)
                if p.is_absolute() or '..' in p.parts or len(p.parts) < 2:
                    raise ValueError('Unsafe binary-kit entry: ' + name)
                relative = str(PurePosixPath(*p.parts[1:]))
                permitted = (relative in {'README.rst', 'VALIDATION.rst', 'AUTHORS.rst',
                    'CITATION.cff', 'CITATIONS.bib', 'THIRD_PARTY_NOTICES.rst',
                    'requirements-validation.txt', 'SHA256SUMS', wheels[0].name}
                    or relative in PUBLIC_HELPERS
                    or p.parts[1] == 'docs' and p.suffix in ('.md', '.rst'))
                if not permitted:
                    raise ValueError('Unexpected binary-kit entry: ' + name)
                if check_rst and p.suffix == '.rst':
                    from docutils.core import publish_doctree
                    publish_doctree(archive.read(name).decode(), source_path='binary-kit/' + relative,
                                   settings_overrides={'halt_level': 2, 'report_level': 2})
                    kit_rst_documents += 1
                if relative == wheels[0].name and hashlib.sha256(archive.read(name)).hexdigest() != wheel_report['sha256']:
                    raise ValueError('Kit contains a different wheel')
        versions.append({'version': manifest['version'], 'wheel_sha256': wheel_report['sha256'],
                         'native_modules': len(wheel_report['native_sha256']), 'protected_source_files': 0})
    if not versions:
        raise ValueError('No verified binary downloads')
    rst_documents = 0
    if check_rst:
        from docutils.core import publish_doctree
        for name in paths:
            if name.endswith('.rst'):
                publish_doctree((root/name).read_text(), source_path=name,
                               settings_overrides={'halt_level': 2, 'report_level': 2})
                rst_documents += 1
    return {'status': 'passed', 'files': len(paths), 'downloads': versions,
            'core_source_in_current_tree': False, 'rst_documents_checked': rst_documents,
            'kit_rst_documents_checked': kit_rst_documents,
            'scope': 'current tree; historical commits/branches are retained'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--json', type=Path)
    parser.add_argument('--check-rst', action='store_true', help='Also render RST documents; requires docutils')
    args = parser.parse_args()
    report = json.dumps(check(args.root, args.check_rst), indent=2) + '\n'
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(report)
    print(report, end='')
