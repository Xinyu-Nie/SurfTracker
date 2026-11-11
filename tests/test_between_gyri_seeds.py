import json
from dataclasses import asdict

import numpy as np
import pytest
from scipy import sparse

from swmtrack import cli
from swmtrack.io import read_surface, write_surface, sha256
from swmtrack.open_preprocessing import METHOD_VERSION
from swmtrack.preprocessing import (PreprocessConfig, intersection_seed_domains,
    preprocess, write_branches)
from test_open_preprocessing import sphere


def test_connected_seeds_keep_all_interface_components_and_cannot_cross_gyrus_or_outside_roi():
    atlas=np.array([1,1,2,1,1,1,2,1,3,2,1,2])
    seg=np.array([1,1,1,0,1,1,1,1,1,1,1,0])
    edges=np.array([[0,1],[1,2],[0,3],[3,4],[5,6],[7,8],[8,9],[10,11]])
    graph=sparse.coo_matrix((np.ones(2*len(edges),bool),
        (np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(12,12)).tocsr()
    raw,connected,boundary,_,details=intersection_seed_domains(atlas,seg,graph,[1],[2])
    np.testing.assert_array_equal(np.flatnonzero(raw),[0,1,2,4,5,6,7,9,10])
    np.testing.assert_array_equal(np.flatnonzero(connected),[0,1,2,5,6,10])
    assert boundary[10] and details['retained_components']==3
    reversed_masks=intersection_seed_domains(atlas,seg,graph,[2],[1])
    for a,b in zip((raw,connected,boundary),reversed_masks[:3]):np.testing.assert_array_equal(a,b)


def test_grouped_rois_ignore_internal_parcel_interfaces_and_unknown_neighbors():
    atlas=np.array([11,12,12,21,11,0]);seg=np.ones(6,int)
    graph=sparse.coo_matrix((np.ones(6,bool),([0,1,2,3,4,5],[1,0,3,2,5,4])),shape=(6,6)).tocsr()
    raw,connected,_,_,_=intersection_seed_domains(atlas,seg,graph,[11,12],[21])
    np.testing.assert_array_equal(raw,[1,1,1,1,1,0])
    np.testing.assert_array_equal(connected,[0,0,1,1,0,0])
    _,empty,_,_,details=intersection_seed_domains(atlas,seg,graph,[11,12],[22])
    assert not empty.any() and details['retained_components']==0
    with pytest.raises(ValueError,match='disjoint'):
        intersection_seed_domains(atlas,seg,graph,[11,12],[12,21])


def test_public_preprocess_preserves_traversable_domain_and_endpoints_when_changing_seeds(tmp_path,monkeypatch):
    p,f=sphere(3);n=len(p)
    atlas=np.where(p[:,1]>=0,17,42)
    # One interface-associated band, plus an isolated sulcal cap inside A.
    seg=(((np.abs(p[:,1])<.25)&(np.abs(p[:,2])<.5))|
         ((p[:,2]>.8)&(p[:,1]>.3))).astype(np.int32)
    source=tmp_path/'surface.vtk'
    write_surface(source,p,f,{'aparc':atlas,'source_vertex_id':np.arange(n)})
    cache=tmp_path/'cache';(cache/'intermediates').mkdir(parents=True)
    old_config=asdict(PreprocessConfig());old_config.pop('between_gyri_seed_policy')
    (cache/'preprocessing.json').write_text(json.dumps({'surface_sha256':sha256(source),
        'atlas_array':'aparc','config':old_config,'method_version':METHOD_VERSION}))
    arrays={k:np.zeros(n) for k in ('feature','depth','flux')}
    arrays.update({k:seg for k in ('initial','cleaned','segmentation','sulcal_skeleton','gyral_skeleton')})
    np.savez_compressed(cache/'segmentation.npz',**arrays)
    write_branches(cache/'intermediates/sulcal_branches.raw',[])
    (cache/'intermediates/flux_pruning.json').write_text('[]')
    def forbidden(*args,**kwargs):raise AssertionError('A skeleton/patch gate was used')
    monkeypatch.setattr('swmtrack.preprocessing.make_patches',forbidden)
    saved={}
    for policy in ('connected-intersection','intersection'):
        config=PreprocessConfig(between_gyri_seed_policy=policy)
        out=tmp_path/policy
        meta=preprocess(source,out,config=config,labels_a=[17],labels_b=[42],reuse_segmentation=cache)
        actual,faces,d=read_surface(out/'preprocessed.world.vtk');saved[policy]=d
        np.testing.assert_array_equal(actual,p);np.testing.assert_array_equal(faces,f)
        np.testing.assert_array_equal(d['source_vertex_id'],np.arange(n))
        np.testing.assert_array_equal(d['sulcal_segmentation'],seg)
        assert meta['tracking_ready'] and meta['seed_policy']==policy
    a,b=saved['connected-intersection'],saved['intersection']
    assert a['seed'].any() and np.any(b['seed']&~a['seed'])
    for field in ('roi','gyrus_a','gyrus_b'):np.testing.assert_array_equal(a[field],b[field])
    assert np.all(a['roi']) and np.all(a['roi'][b['seed'].astype(bool)])
    np.testing.assert_array_equal(b['seed'],seg==1)
    np.testing.assert_array_equal(a['gyrus_a'],(seg==0)&(atlas==17))
    np.testing.assert_array_equal(a['gyrus_b'],(seg==0)&(atlas==42))


@pytest.mark.parametrize('policy',['connected-intersection','intersection','legacy-patches'])
def test_cli_exposes_between_gyri_seed_policy_and_defaults(tmp_path,monkeypatch,policy):
    captured={}
    def capture(*args,**kwargs):captured.update(kwargs);return {}
    monkeypatch.setattr(cli,'preprocess',capture)
    argv=['preprocess','--surface','unused','--output-dir',str(tmp_path)]
    cli.main(argv)
    assert captured['config'].between_gyri_seed_policy=='connected-intersection'
    cli.main(argv+['--between-gyri-seed-policy',policy])
    assert captured['config'].between_gyri_seed_policy==policy
    captured['config'].validate()


def test_unknown_seed_policy_is_rejected():
    with pytest.raises(ValueError,match='between_gyri_seed_policy'):
        PreprocessConfig(between_gyri_seed_policy='wrong').validate()
