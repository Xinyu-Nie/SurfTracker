"""Check SurfTracker wheel/sdist contents without installing or extracting them.

Only the public package, portable tests, documentation and helper scripts
belong in these archives. A project license is reported separately until the
maintainers select the approved terms; --require-license enforces that gate.
"""
import argparse
import ast
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile
import zipfile


FORBIDDEN_PARTS = {"_native", "__pycache__", ".git", ".venv", "results", "research", "research_tests", "build", "dist",
                   "Surface_tracking_c++", "Surface_tracking_matlab"}
FORBIDDEN_SUFFIXES = (".pyc", ".nbc", ".nbi", ".so", ".dll", ".exe", ".dylib", ".o", ".a",
                      ".nii", ".nii.gz", ".npz", ".npy", ".vtk", ".tck", ".trk", ".mat", ".raw")


def check_archive(path, *, require_license=False):
    path = Path(path)
    wheel = path.name.endswith(".whl")
    if wheel:
        with zipfile.ZipFile(path) as archive:
            files = {n: archive.read(n) for n in archive.namelist() if not n.endswith("/")}
        metadata_paths = [n for n in files if n.endswith(".dist-info/METADATA")]
    else:
        with tarfile.open(path) as archive:
            members = [m for m in archive.getmembers() if m.isfile()]
            if any(m.issym() or m.islnk() for m in archive.getmembers()):
                raise ValueError("Source archive must not contain links")
            files = {m.name: archive.extractfile(m).read() for m in members}
        roots = {PurePosixPath(n).parts[0] for n in files}
        if len(roots) != 1:
            raise ValueError("Source archive must have one project root")
        files = {str(PurePosixPath(n).relative_to(next(iter(roots)))): data for n, data in files.items()}
        metadata_paths = [n for n in files if n == "PKG-INFO"]
    if len(metadata_paths) != 1:
        raise ValueError("Archive must contain one primary package metadata file")
    metadata = BytesParser().parsebytes(files[metadata_paths[0]])
    if str(metadata["Name"]).lower() != "surftracker":
        raise ValueError("Unexpected distribution name")
    for name, data in files.items():
        parts = PurePosixPath(name).parts
        if PurePosixPath(name).is_absolute() or ".." in parts or FORBIDDEN_PARTS.intersection(parts):
            raise ValueError(f"Unexpected/private archive path: {name}")
        if name.lower().endswith(FORBIDDEN_SUFFIXES):
            raise ValueError(f"Binary, data or cache payload is not permitted: {name}")
        if wheel and parts[0] != "swmtrack" and not parts[0].endswith(".dist-info"):
            raise ValueError(f"Unexpected wheel payload: {name}")
        if not wheel and parts[0] not in {"swmtrack", "tests", "scripts", "docs", "examples", ".github"} and parts[0].lower() != "surftracker.egg-info":
            if len(parts) != 1 or not (name in {"README.rst", "README.md", "CHANGELOG.rst", "RELEASE_PREPARATION.rst", "THIRD_PARTY_NOTICES.rst",
                    "pyproject.toml", "MANIFEST.in", ".gitignore", "PKG-INFO", "setup.cfg", "setup.py"}
                    or name.upper().startswith(("LICENSE", "COPYING", "NOTICE", "AUTHORS", "CITATION"))):
                raise ValueError(f"Unexpected source payload: {name}")
        if not wheel and name.startswith("tests/") and name.endswith(".py"):
            for node in ast.walk(ast.parse(data, filename=name)):
                modules = ([node.module or ""] if isinstance(node, ast.ImportFrom) else
                           [x.name for x in node.names] if isinstance(node, ast.Import) else [])
                if any(m == "scripts" or m.startswith("scripts.") for m in modules):
                    raise ValueError(f"Portable test imports a research/script namespace: {name}")
    sources = {n: hashlib.sha256(data).hexdigest() for n, data in files.items()
               if n.startswith("swmtrack/") and n.endswith(".py")}
    if "swmtrack/__init__.py" not in sources:
        raise ValueError("Package source is missing")
    module = ast.parse(files["swmtrack/__init__.py"])
    version = next(ast.literal_eval(node.value) for node in module.body if isinstance(node, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == "__version__" for t in node.targets))
    if version != metadata["Version"]:
        raise ValueError("Source and distribution versions differ")
    license_files = [n for n in files if PurePosixPath(n).name.upper().startswith(("LICENSE", "COPYING"))]
    license_declared = bool(metadata["License-Expression"] or metadata["License"])
    if require_license and not (license_files and license_declared):
        raise ValueError("Approved project license text/metadata is still missing")
    if not wheel:
        required = {"tests/test_recommended_defaults.py", "tests/test_whole_hemisphere.py",
                    "tests/test_portable_helpers.py", "scripts/export_qc.py", "scripts/check_distributions.py",
                    "scripts/verify_artifacts.py", "scripts/synthetic_demo.py", ".github/workflows/ci.yml"}
        if required - files.keys():
            raise ValueError(f"Source archive omits portable files: {sorted(required - files.keys())}")
    return {"path": str(path.resolve()), "distribution_name": metadata["Name"], "version": version, "files": len(files),
            "python_sources": len(sources), "source_sha256": sources,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size,
            "project_license_ready": bool(license_files and license_declared)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path, nargs="+")
    parser.add_argument("--require-license", action="store_true")
    parser.add_argument("--json", type=Path, help="Optional report file")
    args = parser.parse_args()
    reports = [check_archive(p, require_license=args.require_license) for p in args.archives]
    baseline = reports[0]
    if any(r["version"] != baseline["version"] or r["source_sha256"] != baseline["source_sha256"] for r in reports):
        raise ValueError("Wheel/source archives differ in version or package source")
    result = {"status": "passed", "archives": reports}
    output = json.dumps(result, indent=2) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
