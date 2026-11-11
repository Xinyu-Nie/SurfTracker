import json
from pathlib import Path
import numpy as np
import nibabel as nib
import pytest
from swmtrack.workflow import resolve_labels, run_pipeline
from swmtrack.tracking import TrackingConfig
from test_freesurfer import make_subject
from swmtrack import workflow


def test_annotation_names_and_ids_resolve_to_same_nonempty_region():
    table=[{'name':'precentral','id':17,'vertices':10},
           {'name':'unknown','id':0,'vertices':4}, {'name':'empty','id':42,'vertices':0}]
    assert resolve_labels(['precentral'],table) == resolve_labels(['17'],table) == [17]
    for request in (['unknown'],['empty'],['absent'],['precentral','17']):
        with pytest.raises(ValueError):resolve_labels(request,table)


def test_workflow_requires_explicit_alignment_basis_and_endpoint_semantics(tmp_path):
    common=dict(subject_dir='unused',fod='unused',output_dir=tmp_path)
    with pytest.raises(ValueError,match='sh_basis'):
        run_pipeline(**common)
    with pytest.raises(ValueError,match='scanner'):
        run_pipeline(**common,sh_basis='legacy-real-even')
    with pytest.raises(ValueError,match='no --gyrus-b'):
        run_pipeline(**common,sh_basis='legacy-real-even',assume_scanner_aligned=True,
                     tracking_config=TrackingConfig(endpoint_policy='within_gyrus'),gyrus_b=['postcentral'])


def test_resume_checks_saved_stage_files_and_run_configuration(tmp_path,monkeypatch):
    subject,fod,_,_=make_subject(tmp_path)
    calls=[]
    def preprocess(surface,output_dir,**kwargs):
        calls.append('preprocess')
        root=Path(output_dir)
        (root/'preprocessed.world.vtk').write_text('test stage payload')
        (root/'preprocessing.json').write_text(json.dumps({'tracking_ready':True,'patch_count':1}))
    class Prepared:
        metadata={'test':True}
    def prepare(*args):
        calls.append('projection'); return Prepared()
    def track(*args):
        calls.append('tracking'); return [],{},{}
    monkeypatch.setattr(workflow,'preprocess',preprocess)
    monkeypatch.setattr(workflow,'prepare',prepare)
    monkeypatch.setattr(workflow,'save_prepared',lambda path,data:Path(path).write_bytes(b'projection'))
    monkeypatch.setattr(workflow,'load_prepared',lambda path:Prepared())
    monkeypatch.setattr(workflow,'track',track)
    monkeypatch.setattr(workflow,'write_tractograms',lambda *args:{})
    monkeypatch.setattr(workflow,'validate_exports',lambda *args:{'streamlines':0})
    options=dict(subject_dir=subject,fod=fod,output_dir=tmp_path/'run',assume_scanner_aligned=True,
                 hemispheres=('lh',),sh_basis='legacy-real-even')
    assert run_pipeline(**options)['status']=='complete'
    assert calls==['preprocess','projection','tracking']
    calls.clear()
    assert run_pipeline(**options,resume=True)['status']=='complete'
    assert not calls
    with pytest.raises(ValueError,match='nonempty'):run_pipeline(**options)
    with pytest.raises(ValueError,match='parameters'):
        run_pipeline(**options,resume=True,shrink_mm=.75)
    (tmp_path/'run/lh/preprocessing/preprocessed.world.vtk').write_text('tampered')
    with pytest.raises(ValueError,match='integrity'):run_pipeline(**options,resume=True)
    assert not calls
