import csv
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
import pytest
from swmtrack import cli, workflow
from swmtrack.io import read_surface,write_surface,sha256,write_streamlines_vtk
from swmtrack.parcellation import cortical_labels,classify_endpoint_parcels,write_endpoint_classification
from swmtrack.preprocessing import preprocess,PreprocessConfig,write_branches
from swmtrack.open_preprocessing import METHOD_VERSION
from swmtrack.tracking import TrackingConfig,track
from test_preprocessing import segmentation_cache
from test_freesurfer import make_subject
from test_geometry_tracking import prepared_grid


def test_cortical_labels_exclude_named_medial_wall_even_with_positive_id():
    table=[dict(id=9,name='unknown',vertices=12),dict(id=10,name='corpuscallosum',vertices=12),
        dict(id=11,name='medial_wall',vertices=12),dict(id=17,name='precentral',vertices=12),
        dict(id=42,name='absent',vertices=0)]
    assert cortical_labels(table)==[17]
    with pytest.raises(ValueError,match='no nonempty'):cortical_labels(table[:3])


def test_endpoint_labels_are_conservative_endian_safe_and_orientation_preserving(tmp_path):
    # Non-gyral vertices must not vote; disagreeing gyral vertices are ambiguous.
    atlas=np.array([14423100,14423100,9221140,9221140,0],dtype='>i4')
    gyrus=np.array([1,1,1,1,0]);faces=np.array([[0,1,4],[2,3,4],[0,2,4]])
    labels,kinds=classify_endpoint_parcels(atlas,gyrus,faces,np.array([[0,1],[1,0],[0,0],[2,1]]))
    np.testing.assert_array_equal(labels,[[14423100,9221140],[9221140,14423100],[14423100,14423100],[0,9221140]])
    np.testing.assert_array_equal(kinds,[2,2,1,3])
    table=[dict(id=14423100,name='precentral',vertices=2),dict(id=9221140,name='superiorparietal',vertices=2)]
    report,fields=write_endpoint_classification(tmp_path,labels,kinds,[10,11,12,13],[7,8,9,10],table)
    assert (report['same_parcel'],report['cross_parcel'],report['ambiguous_parcel'])==(1,2,1)
    counts=list(csv.DictReader((tmp_path/'parcel_pair_counts.csv').open()))
    assert sum(int(r['tracks']) for r in counts)==4
    assert next(r for r in counts if r['parcel_low_id']=='9221140' and r['parcel_high_id']=='14423100')['tracks']=='2'
    rows=list(csv.DictReader((tmp_path/'endpoint_parcels.csv').open()))
    assert rows[0]['endpoint_a_name']=='precentral' and rows[1]['endpoint_a_name']=='superiorparietal'
    np.testing.assert_array_equal(fields['accepted_trial_id'],[7,8,9,10])
    labels,kinds=classify_endpoint_parcels(atlas,gyrus,faces,np.empty((0,2),int))
    assert labels.shape==(0,2) and len(kinds)==0
    empty=tmp_path/'empty';report,_=write_endpoint_classification(empty,labels,kinds,[],[],table)
    assert report['tracks']==0 and len(list(csv.DictReader((empty/'endpoint_parcels.csv').open())))==0
    with pytest.raises(ValueError,match='invalid'):classify_endpoint_parcels(atlas,gyrus,faces,np.array([[-1,0]]))


def test_whole_masks_respect_cortex_and_preserve_original_annotation(tmp_path):
    surface,cache,_=segmentation_cache(tmp_path)
    p,f,d=read_surface(surface);d['cortex']=np.array([0,1,0,1,1,0]);write_surface(surface,p,f,d)
    meta=json.loads((cache/'preprocessing.json').read_text());meta['surface_sha256']=sha256(surface)
    (cache/'preprocessing.json').write_text(json.dumps(meta))
    meta=preprocess(surface,tmp_path/'whole',labels_a=[17,42],labels_b=None,whole_hemisphere=True,
        config=PreprocessConfig(),reuse_segmentation=cache)
    _,_,out=read_surface(tmp_path/'whole/preprocessed.world.vtk')
    np.testing.assert_array_equal(out['seed'],[0,1,0,0,0,0])
    np.testing.assert_array_equal(out['roi'],[0,1,0,1,1,0])
    np.testing.assert_array_equal(out['gyrus'],[0,0,0,1,1,0])
    np.testing.assert_array_equal(out['aparc'],d['aparc'])
    assert meta['tracking_scope']=='whole-hemisphere' and meta['cortex_mask_used']


def test_merged_hemisphere_tracks_require_both_gyral_endpoints_and_length_limits():
    p=prepared_grid(12);p.gyrus=(p.targets>0).astype(np.int64)
    config=TrackingConfig(endpoint_policy='within_gyrus',n_seeds=200,threads=1,min_length_mm=8,max_length_mm=14)
    paths,report,diag=track(p,config)
    assert len(paths)>0
    ids=diag['accepted_trial_ids']
    assert np.all(diag['endpoint_labels'][ids]==1)
    assert np.all((diag['lengths_mm'][ids]>=8)&(diag['lengths_mm'][ids]<=14))
    assert np.all(diag['stops'][ids]==1)


def test_binary_vtk_cell_scalars_and_polylines_roundtrip_independently(tmp_path):
    paths=[np.array([[1,2,3],[4,5,6]],np.float32),np.array([[7,8,9],[10,11,12],[13,14,15]],np.float32)]
    dest=tmp_path/'tracks.vtk';write_streamlines_vtk(dest,paths,{'parcel':np.array([14423100,9221140]),'length_mm':np.array([3.5,7.25])})
    with dest.open('rb') as f:
        for _ in range(4):f.readline()
        assert f.readline()==b'POINTS 5 float\n'
        np.testing.assert_array_equal(np.frombuffer(f.read(60),'>f4').reshape(-1,3),np.concatenate(paths))
        assert f.readline()==b'\n' and f.readline()==b'LINES 2 7\n'
        np.testing.assert_array_equal(np.frombuffer(f.read(28),'>i4'),[2,0,1,3,2,3,4])
        assert f.readline()==b'\n' and f.readline()==b'CELL_DATA 2\n'
        assert f.readline()==b'SCALARS parcel int 1\n' and f.readline()==b'LOOKUP_TABLE default\n'
        np.testing.assert_array_equal(np.frombuffer(f.read(8),'>i4'),[14423100,9221140])
        assert f.readline()==b'\n' and f.readline()==b'SCALARS length_mm double 1\n'
        assert f.readline()==b'LOOKUP_TABLE default\n'
        np.testing.assert_array_equal(np.frombuffer(f.read(16),'>f8'),[3.5,7.25])
    with pytest.raises(ValueError):write_streamlines_vtk(dest,paths,{'bad':np.array([1])})


def test_whole_public_cli_through_import_masks_projection_classification_and_exports(tmp_path,monkeypatch):
    subject,fod,_,_=make_subject(tmp_path)
    def cached_preprocessing(surface,output_dir,**kwargs):
        points,faces,data=read_surface(surface);n=len(points)
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
    args=['run','--subject-dir',str(subject),'--fod',str(fod),'--output-dir',str(output),'--hemisphere','lh',
        '--sh-basis','legacy-real-even','--assume-scanner-aligned','--shrink-mm','0','--n-seeds','8','--whole-hemisphere']
    cli.main(args)
    report=json.loads((output/'run.json').read_text());assert report['status']=='complete'
    assert report['signature']['whole_hemisphere'] and report['signature']['gyrus_a'] is None
    assert report['hemispheres']['lh']['labels_a']==[197121,394500]
    tracks=output/'lh/tractography'
    assert (tracks/'tracks.vtk').exists() and (tracks/'endpoint_parcels.csv').exists()
    assert report['hemispheres']['lh']['tracking']['endpoint_classification']['tracks']==0
    cli.main(args+['--resume'])
    with pytest.raises(SystemExit):cli.main(args+['--gyrus-a','precentral'])
    with pytest.raises(ValueError,match='stopping'):cli.main(args+['--no-stop-at-targets'])
    with pytest.raises(ValueError,match='within_gyrus'):cli.main(args+['--endpoint-policy','between_gyri'])
