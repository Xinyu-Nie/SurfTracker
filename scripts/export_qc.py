"""Export full-resolution label VTKs from a completed swmtrack run.

This helper only copies geometry/labels. It never reruns segmentation,
registration or tracking, and never changes the pipeline output directory.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from nibabel.freesurfer.io import read_geometry
from swmtrack.io import read_surface, write_surface, sha256


def export_qc(run_dir, output_dir, *, hemisphere="both", inflated=False, subject_dir=None):
    run_dir, output_dir = Path(run_dir).resolve(), Path(output_dir).resolve()
    run = json.loads((run_dir / "run.json").read_text())
    if run.get("status") != "complete":
        raise ValueError("QC export requires a completed run.json")
    hemispheres = list(run["hemispheres"]) if hemisphere == "both" else [hemisphere]
    if not hemispheres or any(h not in run["hemispheres"] for h in hemispheres):
        raise ValueError("Requested hemisphere is absent from the completed run")
    if output_dir == run_dir or run_dir in output_dir.parents:
        raise ValueError("Choose a QC output directory outside the pipeline directory")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("QC output directory is nonempty; choose a new directory")
    annotation = run["signature"]["annotation"]
    imported = json.loads((run_dir / "surfaces/freesurfer_import.json").read_text())
    inputs = []
    for h in hemispheres:
        source = run_dir / h / "preprocessing/preprocessed.world.vtk"
        points, faces, data = read_surface(source)
        roi = np.asarray(data["roi"]).ravel() > 0
        segmentation = np.asarray(data["sulcal_segmentation"]).ravel()
        if not np.isin(segmentation[roi], [0, 1]).all():
            raise ValueError("Tracking domain contains an invalid segmentation class")
        segmentation = np.where(roi, segmentation, 2).astype(np.uint8)
        endpoints = (data["gyrus_a"] > 0).astype(np.uint8) + 2 * (data["gyrus_b"] > 0).astype(np.uint8)
        colors = np.full((len(points), 3), 190, dtype=np.uint8)
        for label in imported["hemispheres"][h]["labels"]:
            colors[data[annotation] == label["id"]] = label["rgba"][:3]
        fields = {**data, "candidate_segmentation": segmentation,
                  "anatomical_roi": data[annotation], "anatomical_roi_RGB": colors,
                  "endpoint_region": endpoints}
        fields["candidate_segmentation_RGB"] = np.array(
            [[230, 179, 89], [89, 140, 217], [180, 180, 180]], dtype=np.uint8)[segmentation]
        fields["seed_RGB"] = np.array([[210, 210, 210], [245, 190, 45]], dtype=np.uint8)[
            (data["seed"] > 0).astype(np.uint8)]
        surfaces = [("white", points)]
        inflated_source = None
        if inflated:
            subject = subject_dir or run["signature"].get("subject_dir")
            if not subject:
                raise ValueError("Inflated export requires --subject-dir or a recorded subject_dir")
            inflated_source = Path(subject) / "surf" / f"{h}.inflated"
            inflated_points, inflated_faces = read_geometry(str(inflated_source))
            # A supplied reflecting RAS transform reverses white-face winding.
            # Face order and vertex identity must still match the native input.
            if inflated_points.shape != points.shape or not np.array_equal(
                    np.sort(inflated_faces, axis=1), np.sort(faces, axis=1)):
                raise ValueError("Inflated and white surfaces have different vertex/face identities")
            surfaces.append(("inflated", inflated_points))
        inputs.append((h, source, faces, fields, surfaces, inflated_source))
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {"run_dir": str(run_dir), "annotation": annotation,
              "segmentation_labels": {"0": "gyrus", "1": "sulcus", "2": "outside tracking domain"},
              "endpoint_labels": {"0": "outside target", "1": "A or merged target", "2": "B"},
              "coordinate_note": "White: FOD-world RAS mm. Inflated: native label QC only; never overlay physical tracks on inflation.",
              "hemispheres": {}}
    for h, source, faces, fields, surfaces, inflated_source in inputs:
        folder = output_dir / h
        folder.mkdir()
        record = {"source_sha256": sha256(source), "files": {}}
        if inflated_source is not None:
            record["inflated_source_sha256"] = sha256(inflated_source)
        for kind, points in surfaces:
            path = folder / f"labels.{kind}.vtk"
            write_surface(path, points, faces, fields)
            read_points, read_faces, read_fields = read_surface(path)
            np.testing.assert_array_equal(read_points, points)
            np.testing.assert_array_equal(read_faces, faces)
            for name, values in fields.items():
                np.testing.assert_array_equal(read_fields[name], values)
            record["files"][kind] = {"path": str(path.relative_to(output_dir)), "sha256": sha256(path),
                                     "roundtrip_exact": True}
        report["hemispheres"][h] = record
    (output_dir / "qc.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hemisphere", choices=["lh", "rh", "both"], default="both")
    parser.add_argument("--inflated", action="store_true", help="Also copy native inflated geometry with the same point labels")
    parser.add_argument("--subject-dir", type=Path, help="Optional override for locating surf/HEMI.inflated")
    args = parser.parse_args()
    print(json.dumps(export_qc(args.run_dir, args.output_dir, hemisphere=args.hemisphere,
                               inflated=args.inflated, subject_dir=args.subject_dir), indent=2))


if __name__ == "__main__":
    main()
