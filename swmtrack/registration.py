"""Explicit physical-coordinate conversion for ITK image-resampling affines."""

from pathlib import Path
import json
import os
import subprocess
import sys
import numpy as np
import nibabel as nib
from .freesurfer import affine_matrix
from .io import sha256


def registration_runtime(python=None):
    """Probe the optional isolated ANTsPy runtime before starting a workflow."""
    interpreter=str(Path(python or sys.executable).absolute())
    result=subprocess.run([interpreter,str(Path(__file__).with_name('_ants_registration.py')),'--probe'],
                          capture_output=True,text=True)
    if result.returncode:
        raise RuntimeError('ANTsPy runtime unavailable. Install swmtrack[registration] or select '
                           '--registration-python. '+result.stderr[-1500:])
    return json.loads(result.stdout)


def validate_diffusion_reference(fod, *, b0=None,dwi=None,bvals=None,mask=None):
    """Validate an explicit b0 or corrected DWI reference without loading all DWI."""
    if (b0 is None)==(dwi is None):raise ValueError('ANTs registration needs exactly one of --b0 or --dwi')
    image=nib.load(b0 or dwi)
    if image.shape[:3]!=fod.shape[:3] or not np.allclose(image.affine,fod.affine,rtol=0,atol=1e-5):
        raise ValueError('Registration reference and FOD must share the same spatial grid and affine')
    if b0 is not None:
        if len(image.shape)!=3 or bvals is not None:raise ValueError('--b0 must be 3D and cannot use --bvals')
        indices=None
    else:
        if len(image.shape)!=4 or bvals is None:raise ValueError('--dwi must be 4D with matching --bvals')
        values=np.loadtxt(bvals).ravel()
        if len(values)!=image.shape[3] or not np.isfinite(values).all() or np.any(values<0):
            raise ValueError('b-values must be finite, nonnegative and match the DWI volume count')
        indices=np.flatnonzero(values<100).tolist()
        if not indices:raise ValueError('No b0 volumes with b-value < 100 s/mm2')
    if mask is not None:
        m=nib.load(mask)
        if m.shape!=image.shape[:3] or not np.allclose(m.affine,image.affine,rtol=0,atol=1e-5):
            raise ValueError('DWI mask must match the registration reference grid')
        data=np.asarray(m.dataobj)
        if not np.isfinite(data).all() or np.any(data<0) or not np.any(data>0):raise ValueError('DWI mask must be finite, nonnegative and nonempty')
    return indices


def run_ants_registration(subject,fod,output,*,b0=None,dwi=None,bvals=None,mask=None,
                          python=None,threads=1,initialization='com'):
    """Rigidly align brain.mgz to b0; move surfaces later, never resample FODs."""
    subject,output=Path(subject),Path(output);output.mkdir(parents=True,exist_ok=True)
    if threads<1 or initialization not in ('com','header'):raise ValueError('Invalid ANTs runtime settings')
    indices=validate_diffusion_reference(nib.load(fod),b0=b0,dwi=dwi,bvals=bvals,mask=mask)
    reference=nib.load(b0 or dwi,keep_file_open=True)
    if indices is None:values=np.asarray(reference.dataobj,dtype=np.float32)
    else:
        values=np.zeros(reference.shape[:3],np.float64)
        for i in indices:values+=np.asarray(reference.dataobj[...,i],dtype=np.float64)
        values=(values/len(indices)).astype(np.float32)
    if not np.isfinite(values).all() or float(values.max())<=float(values.min()):
        raise ValueError('Registration b0 must be finite and nonconstant')
    nib.save(nib.Nifti1Image(values,reference.affine),output/'mean_b0.nii.gz')
    brain=nib.load(subject/'mri/brain.mgz');t1=np.asarray(brain.dataobj,dtype=np.float32)
    if len(brain.shape)!=3 or not np.isfinite(t1).all() or not np.any(t1>0):raise ValueError('brain.mgz must be a finite nonempty 3D brain image')
    nib.save(nib.Nifti1Image(t1,brain.affine),output/'t1_brain.nii.gz')
    if mask is not None:
        m=nib.load(mask);nib.save(nib.Nifti1Image((np.asarray(m.dataobj)>0).astype(np.uint8),m.affine),output/'fixed_mask.nii.gz')
    request={'output':str(output.resolve()),'threads':threads,'initialization':initialization,
             'provided_mask':mask is not None,'b0_indices':indices,'b0_threshold_s_mm2':100,
             'random_seed':20261108,'transform':'Rigid'}
    request_path=output/'request.json';request_path.write_text(json.dumps(request,indent=2)+'\n')
    env=dict(os.environ,ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=str(threads))
    with (output/'ants.log').open('w') as log:
        result=subprocess.run([str(Path(python or sys.executable).absolute()),
            str(Path(__file__).with_name('_ants_registration.py')),'--request',str(request_path)],
            env=env,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'ANTs registration failed; inspect {output/"ants.log"}')
    report=json.loads((output/'registration.json').read_text())
    matrix=affine_matrix(np.loadtxt(output/'scanner_ras_to_fod_world.txt'))
    if not np.allclose(matrix[:3,:3].T@matrix[:3,:3],np.eye(3),atol=1e-5) or np.linalg.det(matrix[:3,:3])<0:
        raise ValueError('ANTs result is not a proper rigid point transform')
    report['output_sha256']={p.name:sha256(p) for p in output.iterdir() if p.is_file() and p.name!='registration.json'}
    (output/'registration.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def read_itk_affine_lps(path):
    """Read one text ITK affine: y = A(x-center) + center + translation.

    Returned matrix acts on LPS millimeters. This function does not guess
    whether an external program uses the transform for points or resampling.
    """
    fields = {}
    transforms = []
    for line in Path(path).read_text().splitlines():
        if ":" in line and not line.lstrip().startswith("#"):
            name, value = line.split(":", 1)
            fields[name.strip()] = value.strip()
            if name.strip() == "Transform":
                transforms.append(value.strip())
    if len(transforms) != 1 or not transforms[0].startswith(("MatrixOffsetTransformBase_", "AffineTransform_")):
        raise ValueError("Expected one 3D ITK affine transform")
    parameters = np.fromstring(fields.get("Parameters", ""), sep=" ")
    center = np.fromstring(fields.get("FixedParameters", ""), sep=" ")
    if parameters.shape != (12,) or center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("Expected 12 affine parameters and a finite 3D center")
    matrix = np.eye(4)
    matrix[:3, :3] = parameters[:9].reshape(3, 3)
    matrix[:3, 3] = parameters[9:]+center-matrix[:3, :3]@center
    return affine_matrix(matrix)


def itk_resampling_affine_to_ras_points(path):
    """Convert fixed->moving LPS image pull-map into moving->fixed RAS points.

    For a legacy ANTs image-resampling affine, inversion is essential when
    moving a surface. The experiment verifies this convention against ANTs
    WarpImageMultiTransform on a known asymmetric landmark image.
    """
    flip = np.diag([-1., -1., 1., 1.])
    return affine_matrix(flip @ np.linalg.inv(read_itk_affine_lps(path)) @ flip)
