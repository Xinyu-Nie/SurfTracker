"""Optional compatibility with separately supplied legacy native kernels."""
from pathlib import Path
import json
import os
import platform
import subprocess
import sys
from .io import sha256

TOOL_NAMES = (
    'MeshLaplaceBeltramiEigen', 'Mesh2DTGrid', 'ImplicitCurvatureEstimation',
    'MeshSegmentation', 'CortexSkeleton', 'SkeletonBranchExtraction',
    'SkeletonPruning', 'GeoDist',
)


def bundle_root():
    return Path(__file__).resolve().parent / '_native' / 'linux_x86_64'


def resolve_tools(overrides=None):
    overrides = overrides or {}
    unknown = set(overrides) - set(TOOL_NAMES)
    if unknown:
        raise ValueError(f'Unknown native kernel names: {sorted(unknown)}')
    bundled = set(TOOL_NAMES) - set(overrides)
    if bundled and (sys.platform != 'linux' or platform.machine().lower() not in ('x86_64', 'amd64')):
        raise RuntimeError('Bundled preprocessing requires Linux x86-64. Supply all eight native tool overrides on other platforms.')
    root = bundle_root()
    return {name: Path(overrides[name]).expanduser().resolve() if name in overrides
            else root / 'bin' / name for name in TOOL_NAMES}


def check_native(overrides=None, *, launch=True):
    """Verify payload hashes, executable permissions, and dynamic-loader startup.

    No-argument invocations print usage; their exit status need not be zero.
    Real kernel invocations additionally check their exit status and outputs.
    """
    overrides = overrides or {}
    tools = resolve_tools(overrides)
    manifest_path=bundle_root()/'manifest.json'
    if len(overrides)<len(TOOL_NAMES) and not manifest_path.is_file():
        raise FileNotFoundError('Legacy executables are not distributed with swmtrack. Use the Python backend or supply all eight paths through --tools-json with --backend legacy-native.')
    manifest = json.loads(manifest_path.read_text()) if len(overrides) < len(TOOL_NAMES) else None
    entries = {}
    for name, path in tools.items():
        if not path.is_file():
            raise FileNotFoundError(f'Missing native kernel {name}: {path}; supply --tools-json with --backend legacy-native')
        digest = sha256(path)
        if name not in overrides and digest != manifest['tools'][name]['sha256']:
            raise RuntimeError(f'Native payload checksum mismatch: {name}; reinstall SurfTracker')
        if not os.access(path, os.X_OK):
            raise PermissionError(f'Native kernel is not executable: {path}; check installation permissions and mount options')
        entry = {'path': str(path), 'sha256': digest, 'bundled': name not in overrides}
        if launch:
            try:
                result = subprocess.run([str(path)], capture_output=True, text=True, timeout=15)
            except (OSError, subprocess.TimeoutExpired) as error:
                raise RuntimeError(f'Cannot launch native kernel {name}: {error}') from error
            output = result.stdout + result.stderr
            if result.returncode < 0 or any(token in output for token in
                    ('error while loading shared libraries', 'version `', 'symbol lookup error')):
                raise RuntimeError(f'Native startup failed for {name}: {output.strip()} (exit {result.returncode})')
            entry['usage_exit_code'] = result.returncode
        entries[name] = entry
    return {'platform': sys.platform, 'machine': platform.machine(), 'tools': entries, 'passed': True}
