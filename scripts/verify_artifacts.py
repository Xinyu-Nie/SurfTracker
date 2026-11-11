"""Test extracted source and installed wheel without the research workspace.

Uses the current interpreter's installed runtime/test/build dependencies. No
network access is needed: the wheel is installed with --no-deps --no-index.
This checks artifact independence, not fresh dependency resolution or platforms.
"""
import argparse
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import tempfile

from check_distributions import check_archive


def verify(wheel, sdist, output_dir):
    wheel, sdist, output_dir = Path(wheel).resolve(), Path(sdist).resolve(), Path(output_dir).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Verification output directory must be new or empty")
    wheel_record, source_record = check_archive(wheel), check_archive(sdist)
    if wheel_record["source_sha256"] != source_record["source_sha256"]:
        raise ValueError("Wheel and source archive contain different runtime sources")
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {"status": "running", "distribution_name": wheel_record["distribution_name"], "version": wheel_record["version"],
              "dependency_environment": "current interpreter; no dependency download; no original workspace on PYTHONPATH",
              "artifacts": [wheel_record, source_record], "steps": []}
    report_path = output_dir / "verification.json"

    def save():
        report_path.write_text(json.dumps(report, indent=2) + "\n")

    def command(name, args, cwd, package_path=None):
        env = dict(os.environ, PYTHONPATH=str(package_path) if package_path else "",
                   PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1",
                   NUMBA_CACHE_DIR=str(work / "numba-cache"), PIP_NO_INDEX="1",
                   PIP_DISABLE_PIP_VERSION_CHECK="1")
        log = output_dir / f"{name}.log"
        with log.open("w") as stream:
            result = subprocess.run([sys.executable, *args], cwd=cwd, env=env,
                                    stdout=stream, stderr=subprocess.STDOUT)
        report["steps"].append({"name": name, "returncode": result.returncode,
                                "log": str(log), "command": [sys.executable, *args]})
        save()
        if result.returncode:
            raise RuntimeError(f"{name} failed; inspect {log}")
        print(name + ": passed", flush=True)

    save()
    try:
        with tempfile.TemporaryDirectory(prefix="swmtrack-artifact-check-") as temporary:
            work = Path(temporary)
            extracted = work / "source"
            extracted.mkdir()
            with tarfile.open(sdist) as archive:
                for member in archive.getmembers():
                    path = PurePosixPath(member.name)
                    if path.is_absolute() or ".." in path.parts or member.issym() or member.islnk():
                        raise ValueError("Unsafe source archive member")
                    destination = extracted.joinpath(*path.parts)
                    if member.isdir():
                        destination.mkdir(parents=True, exist_ok=True)
                    elif member.isfile():
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        destination.write_bytes(archive.extractfile(member).read())
                    else:
                        raise ValueError("Unsupported source archive member")
            source = next(extracted.iterdir())
            outside = work / "outside"
            outside.mkdir()
            target = work / "installed"
            source_test = (
                "import pathlib,swmtrack,pytest; "
                f"assert pathlib.Path(swmtrack.__file__).resolve()==pathlib.Path({str(source / 'swmtrack/__init__.py')!r}).resolve(); "
                "print('PACKAGE_UNDER_TEST',swmtrack.__file__); "
                "raise SystemExit(pytest.main(['-q','-p','no:cacheprovider','tests']))")
            command("sdist_tests", ["-c", source_test], source)
            command("install_wheel", ["-m", "pip", "install", "--no-deps", "--no-index", "--target", str(target), str(wheel)], outside)
            wheel_test = (
                "import pathlib,swmtrack,pytest; "
                f"assert pathlib.Path(swmtrack.__file__).resolve()==pathlib.Path({str(target / 'swmtrack/__init__.py')!r}).resolve(); "
                "print('PACKAGE_UNDER_TEST',swmtrack.__file__); "
                f"raise SystemExit(pytest.main(['-q','-p','no:cacheprovider',{str(source / 'tests')!r}]))")
            command("wheel_tests", ["-c", wheel_test], outside, target)
            entry_point_test = (
                "from importlib.metadata import distribution; import pathlib,swmtrack,swmtrack.cli; "
                "d=distribution('SurfTracker'); "
                "assert d.metadata['Name']=='SurfTracker'; assert d.version==swmtrack.__version__; "
                f"assert pathlib.Path(d.locate_file('')).resolve()==pathlib.Path({str(target)!r}).resolve(); "
                "ep=[e for e in d.entry_points if e.group=='console_scripts' and e.name=='swmtrack']; "
                "assert len(ep)==1 and ep[0].value=='swmtrack.cli:main'; "
                "assert ep[0].load() is swmtrack.cli.main; "
                "print('DISTRIBUTION',d.metadata['Name'],d.version,'COMMAND',ep[0].name); "
                "ep[0].load()(['doctor'])")
            command("installed_entry_point", ["-c", entry_point_test], outside, target)
            command("installed_cli", ["-m", "swmtrack", "doctor"], outside, target)
            rebuilt = work / "rebuilt"
            rebuilt.mkdir()
            command("rebuild_from_sdist", ["-c", f"from setuptools.build_meta import build_wheel; build_wheel({str(rebuilt)!r})"], source)
            rebuilt_record = check_archive(next(rebuilt.glob("*.whl")))
            if rebuilt_record["source_sha256"] != wheel_record["source_sha256"]:
                raise ValueError("Rebuilt wheel differs in runtime source")
            report["rebuilt_wheel_source_exact"] = True
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        save()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    verify(args.wheel, args.sdist, args.output_dir)


if __name__ == "__main__":
    main()
