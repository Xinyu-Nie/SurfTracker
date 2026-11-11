import json
from dataclasses import asdict
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from swmtrack import cli, workflow
from swmtrack.io import read_surface, sha256
from swmtrack.open_preprocessing import METHOD_VERSION
from swmtrack.preprocessing import PreprocessConfig, preprocess, write_branches
from test_freesurfer import make_subject
from test_preprocessing import segmentation_cache


def test_within_default_preserves_full_merged_region_and_original_ids_without_patch_gates(tmp_path,monkeypatch):
    surface,cache,stored=segmentation_cache(tmp_path)
    parent=json.loads((cache/'preprocessing.json').read_text())
    parent['config'].pop('within_gyrus_seed_policy')
    (cache/'preprocessing.json').write_text(json.dumps(parent))
    def forbidden(*args,**kwargs):raise AssertionError('A patch gate was used')
    monkeypatch.setattr('swmtrack.preprocessing.make_patches',forbidden)
    meta=preprocess(surface,tmp_path/'merged',labels_a=[17,42],labels_b=None,reuse_segmentation=cache,config=PreprocessConfig())
    p,f,d=read_surface(tmp_path/'merged/preprocessed.world.vtk')
    source,faces,original=read_surface(surface)
    np.testing.assert_array_equal(p,source);np.testing.assert_array_equal(f,faces)
    np.testing.assert_array_equal(d['aparc'],original['aparc'])
    np.testing.assert_array_equal(d['sulcal_segmentation'],stored['segmentation'])
    np.testing.assert_array_equal(d['roi'],[0,1,1,1,1,0])
    np.testing.assert_array_equal(d['seed'],[0,1,1,0,0,0])
    np.testing.assert_array_equal(d['gyrus'],[0,0,0,1,1,0])
    np.testing.assert_array_equal(d['gyrus_a'],d['gyrus'])
    assert not d['gyrus_b'].any() and 'seed_connected' not in d
    assert meta['tracking_ready'] and meta['seed_policy']=='intersection'
    assert meta['labels_a']==[17,42] and meta['labels_b'] is None
    assert meta['endpoint_policy']=='within_gyrus' and meta['native_calls']==[]


def test_legacy_within_patches_remain_an_explicit_option_for_a_merged_region(tmp_path,monkeypatch):
    surface,cache,stored=segmentation_cache(tmp_path);seen={}
    def patches(points,faces,atlas,segmentation,branches,graph,ring2,backend,config,a,b):
        seen.update(a=a,b=b);return [segmentation==1],[{'accepted':True}]
    monkeypatch.setattr('swmtrack.preprocessing.make_patches',patches)
    meta=preprocess(surface,tmp_path/'legacy',labels_a=[17,42],labels_b=None,
        reuse_segmentation=cache,config=PreprocessConfig(within_gyrus_seed_policy='legacy-patches'))
    assert seen=={'a':[17,42],'b':None}
    assert meta['seed_policy']=='legacy-patches' and meta['tracking_ready']


@pytest.mark.parametrize('labels',[[],[0],[99],[17,17]])
def test_invalid_numeric_roi_lists_are_rejected(tmp_path,labels):
    surface,cache,_=segmentation_cache(tmp_path)
    with pytest.raises(ValueError,match='label'):
        preprocess(surface,tmp_path/'invalid',labels_a=labels,labels_b=None,reuse_segmentation=cache)


def test_invalid_within_seed_policy_is_rejected():
    with pytest.raises(ValueError,match='within_gyrus_seed_policy'):
        PreprocessConfig(within_gyrus_seed_policy='connected-intersection').validate()


def test_names_integer_ids_and_mixed_lists_resolve_to_the_same_union():
    table=[{'name':'superiorparietal','id':17,'vertices':3},
           {'name':'inferiorparietal','id':42,'vertices':4}]
    for values in (['superiorparietal','inferiorparietal'],[17,42],['17','42'],['superiorparietal','42']):
        assert workflow.resolve_labels(values,table)==[17,42]
    with pytest.raises(ValueError,match='distinct'):
        workflow.resolve_labels(['superiorparietal',17],table)


@pytest.mark.parametrize('requested',[['precentral','postcentral'],['197121','394500'],['precentral','394500']])
def test_run_cli_merges_names_or_ids_through_real_preparation_and_tracking(tmp_path,monkeypatch,requested):
    subject,fod,_,_=make_subject(tmp_path)
    captured={}
    # Only the synthetic segmentation is cached. Import, public preprocessing,
    # mask construction, projection, tracking and export checks execute normally.
    def cached_preprocessing(surface,output_dir,**kwargs):
        points,faces,data=read_surface(surface);n=len(points)
        captured.update(kwargs)
        cache=tmp_path/'cached_segmentation';(cache/'intermediates').mkdir(parents=True)
        (cache/'preprocessing.json').write_text(json.dumps({'surface_sha256':sha256(surface),
            'atlas_array':'aparc','config':asdict(kwargs['config']),'method_version':METHOD_VERSION}))
        seg=np.array([1,0,1,0],np.int32)
        arrays={key:np.zeros(n) for key in ('feature','depth','flux')}
        arrays.update({key:seg for key in ('initial','cleaned','segmentation','sulcal_skeleton','gyral_skeleton')})
        np.savez_compressed(cache/'segmentation.npz',**arrays)
        write_branches(cache/'intermediates/sulcal_branches.raw',[])
        (cache/'intermediates/flux_pruning.json').write_text('[]')
        return preprocess(surface,output_dir,reuse_segmentation=cache,**kwargs)
    monkeypatch.setattr(workflow,'preprocess',cached_preprocessing)
    output=tmp_path/'run'
    cli.main(['run','--subject-dir',str(subject),'--fod',str(fod),'--output-dir',str(output),
        '--hemisphere','lh','--sh-basis','legacy-real-even','--assume-scanner-aligned',
        '--shrink-mm','0','--n-seeds','8','--gyrus-a',*requested,'--endpoint-policy','within_gyrus'])
    assert captured['labels_a']==[197121,394500] and captured['labels_b'] is None
    assert captured['config'].within_gyrus_seed_policy=='intersection'
    report=json.loads((output/'run.json').read_text());assert report['status']=='complete'
    _,_,d=read_surface(output/'lh/preprocessing/preprocessed.world.vtk')
    np.testing.assert_array_equal(d['aparc'],[197121,394500,197121,394500])
    np.testing.assert_array_equal(d['roi'],np.ones(4))
    np.testing.assert_array_equal(d['seed'],[1,0,1,0])
    np.testing.assert_array_equal(d['gyrus'],[0,1,0,1])
    assert (output/'lh/tractography/tracks.tck').is_file()
    assert report['hemispheres']['lh']['tracking']['config']['endpoint_policy']=='within_gyrus'


def test_preprocess_cli_accepts_numeric_lists_and_explicit_policy(tmp_path,monkeypatch):
    seen={}
    def capture(*args,**kwargs):seen.update(kwargs);return {}
    monkeypatch.setattr(cli,'preprocess',capture)
    base=['preprocess','--surface','unused','--output-dir',str(tmp_path),
          '--endpoint-policy','within_gyrus','--gyrus-a','9221140','14433500']
    cli.main(base)
    assert seen['labels_a']==[9221140,14433500] and seen['labels_b'] is None
    assert seen['config'].within_gyrus_seed_policy=='intersection'
    cli.main(base+['--within-gyrus-seed-policy','legacy-patches'])
    assert seen['config'].within_gyrus_seed_policy=='legacy-patches'
