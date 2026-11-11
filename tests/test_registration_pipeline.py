import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest

from swmtrack import registration,workflow
from swmtrack.io import sha256,write_tractograms
from swmtrack.preprocessing import PreprocessConfig
from test_freesurfer import make_subject


def diffusion(tmp_path):
    a=np.arange(8**3,dtype=np.float32).reshape(8,8,8)+1
    dwi=tmp_path/'dwi.nii';bvals=tmp_path/'bvals';fod=tmp_path/'fod.nii'
    nib.save(nib.Nifti1Image(np.stack((a,7*a,3*a),axis=3),np.eye(4)),dwi)
    nib.save(nib.Nifti1Image(a[...,None],np.eye(4)),fod)
    bvals.write_text('0 1000 50\n')
    return dwi,bvals,fod,a


def test_mean_b0_uses_only_matching_low_b_volumes_without_touching_fod(tmp_path,monkeypatch):
    dwi,bvals,fod,a=diffusion(tmp_path);before=sha256(fod)
    subject=tmp_path/'subject';(subject/'mri').mkdir(parents=True)
    nib.save(nib.MGHImage(a,np.eye(4)),subject/'mri/brain.mgz')
    def fake_worker(command,**kwargs):
        request=json.load(open(command[-1]));out=Path(request['output'])
        np.savetxt(out/'scanner_ras_to_fod_world.txt',np.eye(4))
        (out/'registration.json').write_text(json.dumps(request))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(registration.subprocess,'run',fake_worker)
    report=registration.run_ants_registration(subject,fod,tmp_path/'fit',dwi=dwi,bvals=bvals)
    np.testing.assert_array_equal(nib.load(tmp_path/'fit/mean_b0.nii.gz').get_fdata(),2*a)
    assert report['b0_indices']==[0,2] and sha256(fod)==before


@pytest.mark.parametrize('case',['no-bvals','wrong-count','no-b0','negative','wrong-grid','both'])
def test_invalid_registration_reference_is_rejected(tmp_path,case):
    dwi,bvals,fod,a=diffusion(tmp_path);kw={'dwi':dwi,'bvals':bvals}
    if case=='no-bvals':kw['bvals']=None
    if case=='wrong-count':bvals.write_text('0 1000')
    if case=='no-b0':bvals.write_text('1000 1000 2000')
    if case=='negative':bvals.write_text('0 -1 50')
    if case=='wrong-grid':nib.save(nib.Nifti1Image(np.stack((a,a,a),axis=3),np.diag([2.,1.,1.,1.])),dwi)
    if case=='both':kw['b0']=dwi
    with pytest.raises(ValueError):registration.validate_diffusion_reference(nib.load(fod),**kw)


def test_registration_runtime_keeps_virtual_environment_executable_path(tmp_path,monkeypatch):
    executable=tmp_path/'env/bin/python';executable.parent.mkdir(parents=True);executable.symlink_to('/usr/bin/python3')
    seen=[]
    def probe(command,**kwargs):
        seen.extend(command);return SimpleNamespace(returncode=0,stdout='{"ants_version":"test"}',stderr='')
    monkeypatch.setattr(registration.subprocess,'run',probe)
    registration.registration_runtime(executable)
    assert seen[0]==str(executable)


def test_registration_stage_is_checkpointed_and_transform_is_reused(tmp_path,monkeypatch):
    subject,fod,_,_=make_subject(tmp_path)
    nib.save(nib.load(subject/'mri/orig.mgz'),subject/'mri/brain.mgz')
    b0=tmp_path/'b0.nii';nib.save(nib.Nifti1Image(np.ones((20,20,20),np.float32),np.eye(4)),b0)
    calls=[]
    monkeypatch.setattr(workflow,'registration_runtime',lambda *a:{'ants_version':'test'})
    def fit(subject,fod,out,**kwargs):
        calls.append('registration');np.savetxt(out/'scanner_ras_to_fod_world.txt',np.eye(4))
        (out/'registration.json').write_text('{}')
    def pre(surface,out,**kwargs):
        calls.append('preprocess');assert kwargs['config']==PreprocessConfig.recommended()
        (out/'preprocessed.world.vtk').write_text('stage')
        (out/'preprocessing.json').write_text('{"tracking_ready":true,"patch_count":1}')
    class Prepared:metadata={}
    monkeypatch.setattr(workflow,'run_ants_registration',fit)
    monkeypatch.setattr(workflow,'preprocess',pre)
    monkeypatch.setattr(workflow,'prepare',lambda *a:Prepared())
    monkeypatch.setattr(workflow,'save_prepared',lambda p,o:Path(p).write_text('projection'))
    monkeypatch.setattr(workflow,'load_prepared',lambda p:Prepared())
    monkeypatch.setattr(workflow,'track',lambda *a:([],{},{}))
    monkeypatch.setattr(workflow,'write_tractograms',lambda *a:{})
    monkeypatch.setattr(workflow,'validate_exports',lambda *a:{})
    args=dict(subject_dir=subject,fod=fod,output_dir=tmp_path/'run',hemispheres=('lh',),
              register_ants=True,b0=b0,sh_basis='legacy-real-even')
    assert workflow.run_pipeline(**args)['status']=='complete';assert calls==['registration','preprocess']
    calls.clear();assert workflow.run_pipeline(**args,resume=True)['status']=='complete';assert calls==[]
    (tmp_path/'run/registration/scanner_ras_to_fod_world.txt').write_text('corrupted')
    with pytest.raises(ValueError,match='integrity'):workflow.run_pipeline(**args,resume=True)


def test_alignment_modes_are_explicit_and_registration_inputs_are_not_ignored(tmp_path):
    base=dict(subject_dir='unused',fod='unused',output_dir=tmp_path,sh_basis='legacy-real-even')
    with pytest.raises(ValueError,match='OR'):workflow.run_pipeline(**base,register_ants=True,assume_scanner_aligned=True)
    with pytest.raises(ValueError,match='require register_ants'):workflow.run_pipeline(**base,assume_scanner_aligned=True,b0='unused')


@pytest.mark.parametrize('paths',[[],[np.array([[1.,2.,3.],[4.,5.,6.]],np.float32)]])
def test_optional_vtk_contains_the_same_physical_points(tmp_path,paths):
    ref=nib.Nifti1Image(np.zeros((10,10,10),np.float32),np.eye(4))
    exports=write_tractograms(tmp_path/'tracks',paths,ref,True)
    with open(exports['vtk'],'rb') as f:
        for _ in range(4):f.readline()
        fields=f.readline().split();assert int(fields[1])==sum(map(len,paths))
        if paths:np.testing.assert_array_equal(np.frombuffer(f.read(24),dtype='>f4').reshape(-1,3),paths[0])
