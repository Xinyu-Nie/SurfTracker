"""FreeSurfer-to-tractograms workflow with hash-validated stage checkpoints."""
from dataclasses import asdict
from pathlib import Path
import json
import time
import numpy as np
import nibabel as nib
from . import __version__
from .freesurfer import import_freesurfer, affine_matrix
from .io import sha256, write_tractograms, read_surface
from .preprocessing import preprocess, PreprocessConfig
from .tracking import TrackingConfig, prepare, save_prepared, load_prepared, track
from .validation import validate_exports
from .registration import registration_runtime,validate_diffusion_reference,run_ants_registration
from .parcellation import cortical_labels, classify_endpoint_parcels, write_endpoint_classification


def resolve_labels(requested, table):
    """Resolve annotation names or packed IDs; exclude absent/unknown labels."""
    by_name = {item['name']: item for item in table}
    by_id = {item['id']: item for item in table if item['vertices'] > 0}
    result = []
    for value in requested:
        item = by_name.get(str(value))
        if item is None:
            try:
                item = by_id.get(int(value))
            except (ValueError, TypeError):
                pass
        if item is None or item['id'] <= 0 or item['vertices'] == 0:
            raise ValueError(f'Gyrus label is absent or unknown: {value}')
        result.append(int(item['id']))
    if len(set(result)) != len(result):
        raise ValueError('Gyrus labels must be distinct')
    return result


def run_pipeline(subject_dir, fod, output_dir, *, scanner_to_world=None,
                 assume_scanner_aligned=False, hemispheres=('lh', 'rh'), annotation='aparc',
                 gyrus_a=None, gyrus_b=None, tracking_config=None,whole_hemisphere=False,
                 preprocess_config=None, shrink_mm=.5, interpolation='vertex',
                 sh_basis=None, sh_to_world=None, resume=False,
                 register_ants=False,b0=None,dwi=None,bvals=None,dwi_mask=None,
                 registration_python=None,registration_threads=1,registration_initialization='com'):
    """Run import, preprocessing, projection, tracking, and independent export checks.

    Alignment is explicit. The supported SH basis must also be declared; channel
    count alone cannot identify it. Resume requires identical inputs/config/code
    and verifies every completed stage's saved files before skipping that stage.
    """
    started = time.perf_counter()
    tracking_config = tracking_config or TrackingConfig(endpoint_policy='within_gyrus' if whole_hemisphere else 'between_gyri')
    preprocess_config = preprocess_config or PreprocessConfig.recommended()
    tracking_config.validate(); preprocess_config.validate()
    if sh_basis != 'legacy-real-even':
        raise ValueError('Declare sh_basis="legacy-real-even"; other SH conventions must be converted first')
    if not np.isfinite(shrink_mm) or shrink_mm < 0 or interpolation not in ('vertex', 'centroid'):
        raise ValueError('Invalid shrink distance or interpolation')
    if sum((scanner_to_world is not None,bool(assume_scanner_aligned),bool(register_ants)))!=1:
        raise ValueError('Supply scanner_to_world, assume_scanner_aligned=True, OR register_ants=True')
    if not register_ants and any(v is not None for v in (b0,dwi,bvals,dwi_mask,registration_python)):
        raise ValueError('DWI/b0 and ANTs runtime inputs require register_ants=True')
    if registration_threads<1 or registration_initialization not in ('com','header'):
        raise ValueError('Invalid registration settings')
    scanner_to_world = affine_matrix(np.eye(4) if scanner_to_world is None else scanner_to_world)
    if sh_to_world is not None:
        sh_to_world = np.asarray(sh_to_world, dtype=float)
        if sh_to_world.shape != (3,3) or not np.isfinite(sh_to_world).all() or not np.allclose(sh_to_world.T@sh_to_world, np.eye(3), atol=1e-8):
            raise ValueError('sh_to_world must be an orthogonal 3 x 3 direction transform')
    if not hemispheres or len(set(hemispheres)) != len(hemispheres) or any(h not in ('lh','rh') for h in hemispheres):
        raise ValueError('Select lh and/or rh once each')
    if not annotation or Path(annotation).name != annotation:
        raise ValueError('annotation must be a FreeSurfer annotation name')
    if tracking_config.endpoint_policy == 'unrestricted':
        raise ValueError('The preprocessing workflow requires between_gyri or within_gyrus ROI selection')
    if whole_hemisphere:
        if gyrus_a is not None or gyrus_b is not None:
            raise ValueError('Whole-hemisphere tracking selects all cortical parcels; omit --gyrus-a/b')
        if tracking_config.endpoint_policy!='within_gyrus' or not tracking_config.validate():
            raise ValueError('Whole-hemisphere tracking requires within_gyrus with stopping at gyral targets')
        if preprocess_config.within_gyrus_seed_policy!='intersection':
            raise ValueError('Whole-hemisphere tracking requires intersection seeds')
    elif tracking_config.endpoint_policy == 'within_gyrus':
        gyrus_a = ('precentral',) if gyrus_a is None else gyrus_a
        if not gyrus_a or gyrus_b is not None:
            raise ValueError('within_gyrus requires one or more --gyrus-a labels and no --gyrus-b')
    else:
        gyrus_a = ('precentral',) if gyrus_a is None else gyrus_a
        gyrus_b = ('postcentral',) if gyrus_b is None else gyrus_b
        if not gyrus_a or not gyrus_b:
            raise ValueError('between_gyri requires two nonempty label groups')
    subject, fod, output = Path(subject_dir).resolve(), Path(fod).resolve(), Path(output_dir).resolve()
    relative = ['mri/orig.mgz']
    for hemi in hemispheres:
        relative += [f'surf/{hemi}.white', f'surf/{hemi}.pial', f'label/{hemi}.{annotation}.annot']
        relative += [name for name in (f'surf/{hemi}.{field}' for field in ('curv','sulc','thickness')) if (subject/name).exists()]
        if (subject/f'label/{hemi}.cortex.label').exists():
            relative.append(f'label/{hemi}.cortex.label')
    fod_image = nib.load(fod)
    from .harmonics import degree_from_count
    if len(fod_image.shape) != 4:
        raise ValueError('FOD must be a 4D SH coefficient image')
    degree_from_count(fod_image.shape[-1])
    registration=None
    if register_ants:
        validate_diffusion_reference(fod_image,b0=b0,dwi=dwi,bvals=bvals,mask=dwi_mask)
        relative.append('mri/brain.mgz')
        registration={'runtime':registration_runtime(registration_python),'threads':registration_threads,
            'initialization':registration_initialization,'method':'Rigid','b0_threshold_s_mm2':100,'random_seed':20261108,
            'inputs':{k:{'path':str(Path(v).resolve()),'sha256':sha256(v)}
                      for k,v in dict(b0=b0,dwi=dwi,bvals=bvals,mask=dwi_mask).items() if v is not None}}
    signature = {
        'version': __version__, 'subject_dir': str(subject), 'fod': str(fod),
        'source_sha256': {p.name: sha256(p) for p in Path(__file__).parent.glob('*.py')},
        'input_sha256': {name: sha256(subject/name) for name in relative},
        'fod_sha256': sha256(fod), 'hemispheres': list(hemispheres), 'annotation': annotation,
        'scanner_to_world': None if register_ants else scanner_to_world.tolist(), 'assume_scanner_aligned': assume_scanner_aligned,
        'registration':registration,
        'sh_basis': sh_basis, 'sh_to_world': sh_to_world.tolist() if sh_to_world is not None else None,
        'gyrus_a': list(gyrus_a) if gyrus_a is not None else None, 'gyrus_b': list(gyrus_b) if gyrus_b is not None else None,
        'whole_hemisphere':bool(whole_hemisphere),
        'tracking': asdict(tracking_config), 'preprocessing': asdict(preprocess_config),
        'shrink_mm': shrink_mm, 'interpolation': interpolation,
    }
    manifest_path = output/'run.json'
    if output.exists() and any(output.iterdir()):
        if not resume or not manifest_path.is_file():
            raise ValueError('Output directory is nonempty; choose a new directory or --resume a matching run')
        previous = json.loads(manifest_path.read_text())
        if previous.get('signature') != signature:
            raise ValueError('Cannot resume: inputs, parameters, or source code changed')
    output.mkdir(parents=True, exist_ok=True)
    checkpoints = output/'.checkpoints'; checkpoints.mkdir(exist_ok=True)
    report = {'signature': signature, 'status': 'running', 'backend': 'python', 'hemispheres': {}}

    def save_report():
        temporary = output/'run.json.tmp'
        temporary.write_text(json.dumps(report,indent=2)+'\n')
        temporary.replace(manifest_path)

    def stage(name, directory, action):
        checkpoint = checkpoints/(name+'.json')
        if resume and checkpoint.exists():
            record = json.loads(checkpoint.read_text())
            if not record['files'] or any(not (output/rel).is_file() or sha256(output/rel) != digest for rel,digest in record['files'].items()):
                raise ValueError(f'Resume checkpoint failed file integrity: {name}')
            print(f'{name}: reused verified outputs', flush=True)
            return
        # A previous interrupted stage has no checkpoint. Clear only that stage's
        # managed directory to prevent stale patch files being reused.
        if directory.exists():
            import shutil
            shutil.rmtree(directory)
        directory.mkdir(parents=True)
        print(f'{name}: running', flush=True)
        action()
        files = {str(path.relative_to(output)): sha256(path) for path in sorted(directory.rglob('*')) if path.is_file()}
        checkpoint.write_text(json.dumps({'files':files},indent=2)+'\n')
        print(f'{name}: complete', flush=True)

    save_report()
    try:
        if register_ants:
            aligned=output/'registration'
            stage('registration',aligned,lambda:run_ants_registration(subject,fod,aligned,b0=b0,dwi=dwi,
                bvals=bvals,mask=dwi_mask,python=registration_python,threads=registration_threads,
                initialization=registration_initialization))
            scanner_to_world=affine_matrix(np.loadtxt(aligned/'scanner_ras_to_fod_world.txt'))
            report['registration']=json.loads((aligned/'registration.json').read_text())
        report['scanner_ras_to_fod_world']=scanner_to_world.tolist()
        imported = output/'surfaces'
        stage('import', imported, lambda: import_freesurfer(subject, fod, imported,
            scanner_to_world=None if assume_scanner_aligned else scanner_to_world,
            assume_scanner_aligned=assume_scanner_aligned, hemispheres=hemispheres, annotation=annotation))
        imported_metadata = json.loads((imported/'freesurfer_import.json').read_text())
        # Resolve every hemisphere before running a potentially expensive stage.
        labels = {}
        for hemi in hemispheres:
            table = imported_metadata['hemispheres'][hemi]['labels']
            a = cortical_labels(table) if whole_hemisphere else resolve_labels(gyrus_a,table)
            b = resolve_labels(gyrus_b,table) if gyrus_b is not None else None
            if b is not None and set(a)&set(b):
                raise ValueError('Between-gyrus groups must be disjoint')
            labels[hemi] = a,b
        for hemi in hemispheres:
            a,b = labels[hemi]
            pre = output/hemi/'preprocessing'
            pre_options={'whole_hemisphere':True} if whole_hemisphere else {}
            stage(hemi+'_preprocess', pre, lambda: preprocess(imported/f'{hemi}.white.world.vtk',pre,
                config=preprocess_config, labels_a=a, labels_b=b, atlas_array=annotation,**pre_options))
            pre_meta = json.loads((pre/'preprocessing.json').read_text())
            if not pre_meta['tracking_ready']:
                raise ValueError(f'{hemi}: no qualifying seeds/endpoints; inspect preprocessing.json and the selected ROIs')
            projection = output/hemi/'projection'
            def project():
                prepared = prepare(pre/'preprocessed.world.vtk',fod,shrink_mm,interpolation,sh_to_world,
                                   projection/'tracking_surface.world.vtk')
                save_prepared(projection/'prepared.npz',prepared)
            stage(hemi+'_projection', projection, project)
            tracks = output/hemi/'tractography'
            def propagate():
                prepared = load_prepared(projection/'prepared.npz')
                streamlines, tracking, diagnostics = track(prepared,tracking_config)
                cell_data=None
                if whole_hemisphere:
                    _,pf,pd=read_surface(projection/'tracking_surface.world.vtk')
                    accepted=diagnostics['accepted_trial_ids']
                    endpoint_parcels,kinds=classify_endpoint_parcels(pd[annotation],pd['gyrus'],pf,
                        diagnostics['endpoint_faces'][accepted])
                    classification,cell_data=write_endpoint_classification(tracks,endpoint_parcels,kinds,
                        diagnostics['lengths_mm'][accepted],accepted,imported_metadata['hemispheres'][hemi]['labels'])
                    tracking['endpoint_classification']=classification
                export_options={'vtk_cell_data':cell_data} if whole_hemisphere else {}
                tracking['exports'] = write_tractograms(tracks/'tracks', streamlines, fod_image, True,**export_options)
                tracking['export_validation'] = validate_exports(tracks/'tracks',streamlines,fod_image,prepared,diagnostics)
                tracking['preparation'] = prepared.metadata
                np.savez_compressed(tracks/'diagnostics.npz',**diagnostics)
                (tracks/'tracking.json').write_text(json.dumps(tracking,indent=2)+'\n')
            stage(hemi+'_tracking', tracks, propagate)
            tracking = json.loads((tracks/'tracking.json').read_text())
            report['hemispheres'][hemi] = {
                'labels_a':a, 'labels_b':b, 'patch_count':pre_meta['patch_count'],
                'tracking': tracking, 'tck':str(tracks/'tracks.tck'), 'trk':str(tracks/'tracks.trk'),'vtk':str(tracks/'tracks.vtk'),
            }
            save_report()
        report['status'] = 'complete'
        report['elapsed_seconds_this_invocation'] = time.perf_counter()-started
        save_report()
        return report
    except Exception as error:
        report.update(status='failed',error=str(error),elapsed_seconds_this_invocation=time.perf_counter()-started)
        save_report()
        raise
