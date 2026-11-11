import numpy as np
import pytest
import json
from dataclasses import asdict
from scipy import sparse
from swmtrack.preprocessing import (write_mni_obj,read_branches,write_branches,adjacency,
    expand_ring,components,flux_prune,regular_subset,reconstruct_regions,validate_full_surface)
from swmtrack.preprocessing import preprocess,PreprocessConfig,connected_join
from swmtrack.io import read_legacy_obj,write_surface,read_surface,sha256
from swmtrack.open_preprocessing import METHOD_VERSION


def octahedron():
    points=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.],[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]])
    faces=np.array([[0,1,2],[0,2,3],[0,3,4],[0,4,1],[5,2,1],[5,3,2],[5,4,3],[5,1,4]])
    return points,faces


def test_mni_interchange_retains_world_coordinates_and_identity(tmp_path):
    points,faces=octahedron();points=points+np.array([11.123456789,-20.,3.])
    path=tmp_path/'surface.obj';write_mni_obj(path,points,faces)
    actual,indices,normals=read_legacy_obj(path)
    np.testing.assert_array_equal(actual,points);np.testing.assert_array_equal(indices,faces)
    np.testing.assert_allclose(np.linalg.norm(normals,axis=1),1)


def test_branch_file_zero_based_roundtrip_and_malformed_input(tmp_path):
    path=tmp_path/'branches.raw';branches=[np.array([0,1,3]),np.array([2,4])]
    write_branches(path,branches)
    for a,b in zip(read_branches(path,5),branches):np.testing.assert_array_equal(a,b)
    np.array([4,0,1],dtype='<i4').tofile(path)
    with pytest.raises(ValueError,match='Malformed'):read_branches(path,5)
    np.array([1,5],dtype='<i4').tofile(path)
    with pytest.raises(ValueError,match='outside'):read_branches(path,5)


def test_two_ring_includes_triangle_neighbors_without_self():
    _,faces=octahedron();graph=adjacency(faces,6);ring2=expand_ring(graph)
    assert graph[0,5]==0 and ring2[0,5]
    assert ring2[0,1] and not np.any(ring2.diagonal())


def test_subset_removes_pinched_and_isolated_regions():
    _,faces=octahedron()
    mask=np.array([1,1,0,1,0,0],dtype=bool)
    assert not regular_subset(mask,faces).any()
    hemisphere=np.array([1,1,1,1,1,0],dtype=bool)
    np.testing.assert_array_equal(regular_subset(hemisphere,faces),hemisphere)


def test_component_reconstruction_flips_regions_without_retained_skeletons():
    graph=sparse.diags([np.ones(5),np.ones(5)],[-1,1],shape=(6,6)).astype(bool).tocsr()
    result=reconstruct_regions(np.array([1,1,0,0,1,1]),[np.array([0])],[np.array([2])],graph)
    np.testing.assert_array_equal(result,[1,1,0,0,0,0])


def test_flux_postpruning_retains_removal_before_threshold_order_and_guards():
    ring=sparse.coo_matrix((np.ones(4,dtype=bool),([1,2,3,4],[2,1,4,3])),shape=(6,6)).tocsr()
    branches=[np.array([0,1]),np.array([2,3]),np.array([4,5])]
    result,history=flux_prune(branches,np.array([-.1,-.1,-20.,-20.,-6.,-6.]),ring,-10,100)
    assert len(result)==1
    np.testing.assert_array_equal(result[0],[2,3])
    assert history[-1]['score']==-12
    result,history=flux_prune([np.array([0,1])],np.ones(6),ring,-10,100)
    assert len(result)==1 and history[-1]['stop']=='no_leaf_guard'


def test_preprocessing_rejects_open_meshes_and_inconsistent_winding():
    points,faces=octahedron();validate_full_surface(points,faces)
    with pytest.raises(ValueError,match='closed'):validate_full_surface(points,faces[:-1])
    bad=faces.copy();bad[0]=bad[0,[0,2,1]]
    with pytest.raises(ValueError,match='winding'):validate_full_surface(points,bad)


def segmentation_cache(tmp_path):
    points,faces=octahedron();surface=tmp_path/'surface.vtk'
    write_surface(surface,points,faces,{'aparc':np.array([0,17,17,42,42,0])})
    cache=tmp_path/'cache';(cache/'intermediates').mkdir(parents=True)
    (cache/'preprocessing.json').write_text(json.dumps({'surface_sha256':sha256(surface),
        'atlas_array':'aparc','config':asdict(PreprocessConfig()),
        'method_version':METHOD_VERSION}))
    arrays={key:np.arange(6,dtype=float)/10 for key in ('feature','depth','flux')}
    arrays.update({key:np.array([0,1,1,0,0,0]) for key in
        ('initial','cleaned','segmentation','sulcal_skeleton','gyral_skeleton')})
    np.savez_compressed(cache/'segmentation.npz',**arrays)
    write_branches(cache/'intermediates/sulcal_branches.raw',[np.array([1,2])])
    (cache/'intermediates/flux_pruning.json').write_text('[]')
    return surface,cache,arrays


def test_segmentation_cache_reuses_results_without_native_kernels(tmp_path):
    surface,cache,arrays=segmentation_cache(tmp_path)
    result=preprocess(surface,tmp_path/'new',reuse_segmentation=cache,stop_after='segmentation',
                       config=PreprocessConfig(patch_include_fraction=.1))
    points,faces,data=read_surface(tmp_path/'new/preprocessed.world.vtk')
    np.testing.assert_array_equal(data['k1'],arrays['feature'])
    np.testing.assert_array_equal(data['sulcal_segmentation'],arrays['segmentation'])
    assert result['native_calls']==[]


def test_segmentation_cache_rejects_mismatched_surface_and_parameters(tmp_path):
    surface,cache,_=segmentation_cache(tmp_path)
    with pytest.raises(ValueError,match='parameters differ'):
        preprocess(surface,tmp_path/'out1',reuse_segmentation=cache,stop_after='segmentation',
                   config=PreprocessConfig(feature_smoothing=.7))
    p,f,d=read_surface(surface);write_surface(surface,p+[1,0,0],f,d)
    with pytest.raises(ValueError,match='different surface'):
        preprocess(surface,tmp_path/'out2',reuse_segmentation=cache,stop_after='segmentation')


@pytest.mark.parametrize('reverse_a,reverse_b',[(False,False),(True,False),(False,True),(True,True)])
def test_connected_join_handles_all_endpoint_orientations(reverse_a,reverse_b):
    ring=sparse.coo_matrix((np.ones(2,dtype=bool),([1,2],[2,1])),shape=(4,4)).tocsr()
    a=np.array([0,1]);b=np.array([2,3])
    joined=connected_join(a[::-1] if reverse_a else a,b[::-1] if reverse_b else b,ring)
    np.testing.assert_array_equal(joined,[0,1,2,3])
    assert connected_join(np.array([0]),np.array([3]),ring) is None
