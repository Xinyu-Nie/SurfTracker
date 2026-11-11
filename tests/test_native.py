import json
from pathlib import Path
import pytest
from swmtrack import native


@pytest.mark.skipif(not (native.bundle_root()/'manifest.json').exists(),reason='Optional local legacy payload is not distributed')
def test_bundled_payload_integrity_and_all_eight_executables_start():
    report = native.check_native()
    assert report['passed'] and set(report['tools']) == set(native.TOOL_NAMES)
    assert all(item['bundled'] and 'usage_exit_code' in item for item in report['tools'].values())


@pytest.mark.skipif(not (native.bundle_root()/'manifest.json').exists(),reason='Optional local legacy payload is not distributed')
def test_corrupt_bundle_is_rejected_before_execution(tmp_path, monkeypatch):
    root = native.bundle_root()
    (tmp_path/'bin').mkdir()
    manifest = json.loads((root/'manifest.json').read_text())
    name = native.TOOL_NAMES[0]
    (tmp_path/manifest['tools'][name]['file']).write_bytes(b'corrupted executable')
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    monkeypatch.setattr(native,'bundle_root',lambda:tmp_path)
    with pytest.raises(RuntimeError,match='checksum'):
        native.check_native()


def test_unknown_tool_override_and_unsupported_platform_rejected(tmp_path, monkeypatch):
    with pytest.raises(ValueError,match='Unknown'):
        native.resolve_tools({'typo':'/tmp/tool'})
    monkeypatch.setattr(native.sys,'platform','darwin')
    with pytest.raises(RuntimeError,match='Linux x86-64'):
        native.resolve_tools()
    custom = {name:str(tmp_path/name) for name in native.TOOL_NAMES}
    assert native.resolve_tools(custom) == {name:Path(path).resolve() for name,path in custom.items()}


def test_binary_free_installation_has_clear_legacy_diagnostic(tmp_path,monkeypatch):
    monkeypatch.setattr(native,'bundle_root',lambda:tmp_path)
    monkeypatch.setattr(native.sys,'platform','linux')
    monkeypatch.setattr(native.platform,'machine',lambda:'x86_64')
    with pytest.raises(FileNotFoundError,match='not distributed'):
        native.check_native()
    from swmtrack.open_preprocessing import check_python
    assert check_python()['passed']
