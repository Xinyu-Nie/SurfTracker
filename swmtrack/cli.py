import argparse
from dataclasses import fields, replace
from pathlib import Path
import json
import numpy as np
import nibabel as nib
from .io import convert_legacy, write_tractograms
from .tracking import TrackingConfig, prepare, save_prepared, load_prepared, track
from .freesurfer import import_freesurfer
from .preprocessing import preprocess, PreprocessConfig
from .native import check_native
from .workflow import run_pipeline


def add_tracking_options(parser, *, workflow=False):
    defaults = TrackingConfig()
    for name in ("seed", "threads", "n_seeds", "samples_per_face", "max_points_per_half", "batch_size"):
        parser.add_argument("--"+name.replace("_", "-"), type=int, default=getattr(defaults, name))
    for name in ("seeds_per_mm2", "step_mm", "max_turn_deg_per_mm", "angular_resolution_deg", "min_fod", "min_length_mm", "max_length_mm"):
        parser.add_argument("--"+name.replace("_", "-"), type=float, default=getattr(defaults, name))
    parser.add_argument("--seed-mode", choices=["area", "centroid"], default="area")
    parser.add_argument("--endpoint-policy", choices=["between_gyri", "within_gyrus"] if workflow else ["between_gyri", "within_gyrus", "unrestricted"], default=None if workflow else "between_gyri")
    parser.add_argument("--stop-at-targets", action=argparse.BooleanOptionalAction, default=None)


def add_preprocessing_options(parser):
    parser.add_argument("--preprocess-profile", choices=["recommended", "default", "agreement-candidate", "previous-defaults"],
                        default="recommended", help="Validated recommended settings by default; previous-defaults restores the earlier cotan profile")
    for field in fields(PreprocessConfig):
        default=getattr(PreprocessConfig(),field.name)
        options={}
        if field.name=="between_gyri_seed_policy":
            options=dict(choices=["connected-intersection","intersection","legacy-patches"],
                         help="Between-gyri seeds; default connected-intersection.")
        elif field.name=="within_gyrus_seed_policy":
            options=dict(choices=["intersection","legacy-patches"],
                         help="Within-gyrus seeds in one or merged ROIs; default intersection.")
        parser.add_argument("--"+field.name.replace("_","-"),type=type(default),default=argparse.SUPPRESS,**options)


def preprocessing_config(args):
    config=PreprocessConfig() if args.preprocess_profile=="previous-defaults" else PreprocessConfig.recommended()
    return replace(config,**{f.name:getattr(args,f.name) for f in fields(PreprocessConfig) if hasattr(args,f.name)})


def main(argv=None):
    parser = argparse.ArgumentParser(description="World-space superficial white matter surface tracking")
    commands = parser.add_subparsers(dest="command", required=True)
    doctor = commands.add_parser("doctor", help="Verify the selected preprocessing runtime")
    doctor.add_argument("--backend", choices=["python","legacy-native"], default="python")
    doctor.add_argument("--tools-json", help="Optional JSON mapping of kernel names to replacement executables")
    run = commands.add_parser("run", help="FreeSurfer + FOD to world-space surfaces, masks, and TCK/TRK")
    run.add_argument("--subject-dir", required=True)
    run.add_argument("--fod", required=True)
    run.add_argument("--output-dir", required=True)
    run.add_argument("--hemisphere", choices=["lh", "rh", "both"], default="both")
    run.add_argument("--annotation", default="aparc")
    run.add_argument("--gyrus-a", nargs="+", help="One or more annotation names or packed IDs; default precentral; within_gyrus merges this group")
    run.add_argument("--gyrus-b", nargs="+", help="Second group for between_gyri; default postcentral")
    run.add_argument("--whole-hemisphere", action="store_true", help="Seed all cortical sulci, terminate both ends in the merged gyral target, and classify endpoint parcels; omit --gyrus-a/b")
    run.add_argument("--sh-basis", choices=["legacy-real-even"], required=True,
                     help="Explicitly declare the supported SH convention; see README")
    run.add_argument("--sh-to-world", help="Orthogonal 3 x 3 SH direction transform")
    run.add_argument("--shrink-mm", type=float, default=.5)
    run.add_argument("--interpolation", choices=["vertex", "centroid"], default="vertex")
    run.add_argument("--tools-json")
    run.add_argument("--resume", action="store_true", help="Reuse hash-verified completed stages with identical inputs/configuration")
    alignment = run.add_mutually_exclusive_group(required=True)
    alignment.add_argument("--scanner-to-world", help="Text 4 x 4 scanner RAS to FOD RAS affine")
    alignment.add_argument("--assume-scanner-aligned", action="store_true")
    alignment.add_argument("--register-ants", action="store_true", help="Optional rigid T1-to-b0 registration using ANTsPy")
    reference=run.add_mutually_exclusive_group()
    reference.add_argument("--b0", help="3D mean b0 in the FOD grid, for optional registration")
    reference.add_argument("--dwi", help="Corrected 4D DWI in the FOD grid; supply --bvals to extract mean b0")
    run.add_argument("--bvals", help="b-values matching --dwi volumes")
    run.add_argument("--dwi-mask", help="Optional binary mask in the b0/FOD grid")
    run.add_argument("--registration-python", help="Python interpreter containing ANTsPy; defaults to this interpreter")
    run.add_argument("--registration-threads", type=int, default=1)
    run.add_argument("--registration-initialization", choices=["com","header"], default="com")
    add_preprocessing_options(run)
    pre = commands.add_parser("preprocess", help="Mesh curvature, sulcal/gyral segmentation, skeletons, and sulcal patches")
    pre.add_argument("--surface", required=True)
    pre.add_argument("--output-dir", required=True)
    pre.add_argument("--atlas-array", default="aparc")
    pre.add_argument("--gyrus-a", type=int, nargs="+", default=[14423100], help="One or more packed annotation IDs; within_gyrus merges this group")
    pre.add_argument("--gyrus-b", type=int, nargs="+", default=[1316060])
    pre.add_argument("--endpoint-policy", choices=["between_gyri","within_gyrus"], default="between_gyri")
    pre.add_argument("--tools-json", help="JSON mapping of legacy kernel names to executable paths")
    pre.add_argument("--input-feature", help="Optional float32 per-vertex feature in the selected backend's convention; see README")
    pre.add_argument("--reuse-segmentation", help="Reuse a matching preprocessing output directory when selecting new gyri")
    pre.add_argument("--stop-after", choices=["features","segmentation","patches"], default="patches")
    add_preprocessing_options(pre)
    fs = commands.add_parser("import-freesurfer", help="Import native FreeSurfer surfaces and labels into world-space VTK")
    fs.add_argument("--subject-dir", required=True)
    fs.add_argument("--fod", required=True)
    fs.add_argument("--output-dir", required=True)
    fs.add_argument("--hemisphere", choices=["lh", "rh", "both"], default="both")
    fs.add_argument("--annotation", default="aparc")
    fs.add_argument("--gyrus-a", help="Optional full annotation region name, e.g. precentral")
    fs.add_argument("--gyrus-b", help="Optional full annotation region name, e.g. postcentral")
    transform = fs.add_mutually_exclusive_group(required=True)
    transform.add_argument("--scanner-to-world", help="Text 4 x 4 affine from scanner RAS mm to FOD world RAS mm")
    transform.add_argument("--assume-scanner-aligned", action="store_true", help="Explicitly use scanner RAS as FOD world RAS; verify anatomical alignment")
    convert = commands.add_parser("convert", help="Convert aligned MNI polygon OBJ to a world-space VTK mesh")
    convert.add_argument("--mesh", required=True); convert.add_argument("--fod", required=True)
    convert.add_argument("--output", required=True)
    convert.add_argument("--input-space", choices=["image-mm", "world", "voxel"], default="image-mm")
    convert.add_argument("--mesh-to-world", help="Optional text file containing an explicit 4 x 4 transform")
    for name in ("roi", "seed", "gyrus-a", "gyrus-b", "gyrus"):
        convert.add_argument("--"+name, help="Legacy int32 per-vertex mask")
    prep = commands.add_parser("prepare", help="Build a reusable geometry and Fourier-coefficient cache")
    prep.add_argument("--surface", required=True); prep.add_argument("--fod", required=True)
    prep.add_argument("--output", required=True); prep.add_argument("--roi-output")
    prep.add_argument("--shrink-mm", type=float, default=0.5)
    prep.add_argument("--interpolation", choices=["vertex", "centroid"], default="vertex")
    prep.add_argument("--sh-to-world", help="Orthogonal 3 x 3 direction transform; default SH axes are world RAS axes")
    prep.add_argument("--gyrus-array", default="gyrus", help="One gyrus mask/atlas point array for within_gyrus; default: gyrus")
    prep.add_argument("--gyrus-label", type=int, help="Select one positive atlas label from --gyrus-array")
    tracking = commands.add_parser("track", help="Track from a prepared cache and write both TCK and TRK")
    tracking.add_argument("--cache", required=True); tracking.add_argument("--output", required=True, help="Output prefix")
    add_tracking_options(tracking)
    add_tracking_options(run, workflow=True)
    args = parser.parse_args(argv)
    if args.command == "doctor":
        if args.backend=="python":
            if args.tools_json:parser.error("--tools-json requires --backend legacy-native")
            from .open_preprocessing import check_python
            report=check_python()
        else:report=check_native(json.loads(Path(args.tools_json).read_text()) if args.tools_json else None)
        print(json.dumps(report,indent=2))
    elif args.command == "run":
        if args.whole_hemisphere and (args.gyrus_a is not None or args.gyrus_b is not None):
            parser.error("--whole-hemisphere cannot be combined with --gyrus-a or --gyrus-b")
        args.endpoint_policy=args.endpoint_policy or ("within_gyrus" if args.whole_hemisphere else "between_gyri")
        metadata = run_pipeline(args.subject_dir, args.fod, args.output_dir,
            scanner_to_world=np.loadtxt(args.scanner_to_world) if args.scanner_to_world else None,
            assume_scanner_aligned=args.assume_scanner_aligned,
            hemispheres=("lh", "rh") if args.hemisphere == "both" else (args.hemisphere,),
            annotation=args.annotation, gyrus_a=args.gyrus_a, gyrus_b=args.gyrus_b,whole_hemisphere=args.whole_hemisphere,
            sh_basis=args.sh_basis, sh_to_world=np.loadtxt(args.sh_to_world) if args.sh_to_world else None,
            shrink_mm=args.shrink_mm, interpolation=args.interpolation, resume=args.resume,
            tracking_config=TrackingConfig(**{f.name:getattr(args,f.name) for f in fields(TrackingConfig)}),
            preprocess_config=preprocessing_config(args),
            register_ants=args.register_ants,b0=args.b0,dwi=args.dwi,bvals=args.bvals,dwi_mask=args.dwi_mask,
            registration_python=args.registration_python,registration_threads=args.registration_threads,
            registration_initialization=args.registration_initialization,
            tool_overrides=json.loads(Path(args.tools_json).read_text()) if args.tools_json else None)
        print(json.dumps({"status":metadata["status"],"outputs":{
            hemi:{key:value for key,value in item.items() if key in ("tck","trk","vtk","patch_count")}
            for hemi,item in metadata["hemispheres"].items()}}, indent=2))
    elif args.command == "preprocess":
        config=preprocessing_config(args)
        metadata=preprocess(args.surface,args.output_dir,config=config,labels_a=args.gyrus_a,
            labels_b=args.gyrus_b if args.endpoint_policy=="between_gyri" else None,
            atlas_array=args.atlas_array,input_feature=args.input_feature,stop_after=args.stop_after,
            reuse_segmentation=args.reuse_segmentation,
            tool_overrides=json.loads(Path(args.tools_json).read_text()) if args.tools_json else None)
        print(json.dumps(metadata,indent=2))
    elif args.command == "import-freesurfer":
        metadata = import_freesurfer(args.subject_dir, args.fod, args.output_dir,
            scanner_to_world=np.loadtxt(args.scanner_to_world) if args.scanner_to_world else None,
            assume_scanner_aligned=args.assume_scanner_aligned,
            hemispheres=("lh", "rh") if args.hemisphere == "both" else (args.hemisphere,),
            annotation=args.annotation, gyrus_a=args.gyrus_a, gyrus_b=args.gyrus_b)
        print(json.dumps(metadata, indent=2))
    elif args.command == "convert":
        masks = {name: getattr(args, name) for name in ("roi", "seed", "gyrus_a", "gyrus_b", "gyrus") if getattr(args, name)}
        metadata = convert_legacy(args.mesh, args.fod, args.output, masks, args.input_space,
                                  np.loadtxt(args.mesh_to_world) if args.mesh_to_world else None)
        print(json.dumps(metadata, indent=2))
    elif args.command == "prepare":
        prepared = prepare(args.surface, args.fod, args.shrink_mm, args.interpolation,
                           np.loadtxt(args.sh_to_world) if args.sh_to_world else None, args.roi_output,
                           gyrus_array=args.gyrus_array, gyrus_label=args.gyrus_label)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        save_prepared(args.output, prepared)
        print(json.dumps(prepared.metadata, indent=2))
    else:
        prepared = load_prepared(args.cache)
        config = TrackingConfig(**{f.name: getattr(args, f.name) for f in fields(TrackingConfig)})
        streamlines, report, diagnostics = track(prepared, config)
        reference = nib.load(prepared.metadata["reference_fod"])
        report["exports"] = write_tractograms(args.output, streamlines, reference)
        report["preparation"] = prepared.metadata
        Path(str(args.output)+".json").write_text(json.dumps(report, indent=2)+"\n")
        np.savez(str(args.output)+".diagnostics.npz", **diagnostics)
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
