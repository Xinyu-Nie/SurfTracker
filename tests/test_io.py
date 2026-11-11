import numpy as np
import nibabel as nib
import pytest
from scipy.spatial.transform import Rotation
from swmtrack.io import sample_sh,write_surface,read_surface,write_tractograms
from swmtrack.validation import validate_exports
from swmtrack.geometry import build_geometry
from swmtrack.tracking import Prepared,prepare,track,TrackingConfig


def test_sampling_respects_oblique_anisotropic_affine_and_scaling(tmp_path):
    ijk=np.indices((5,6,7)).astype(float)
    values=(ijk[0]+2*ijk[1]+3*ijk[2]).astype(np.int16)
    affine=np.eye(4);affine[:3,:3]=Rotation.from_euler('xyz',[20,10,30],degrees=True).as_matrix()@np.diag([1.,2.,3.]);affine[:3,3]=[15,-20,9]
    image=nib.Nifti1Image(values[...,None],affine)
    image.header.set_slope_inter(2.,3.)
    path=tmp_path/'scaled.nii';nib.save(image,path);image=nib.load(path)
    voxel=np.array([[1.2,2.3,3.4],[4.,5.,6.]])
    sampled=sample_sh(image,nib.affines.apply_affine(affine,voxel))[:,0]
    np.testing.assert_allclose(sampled,2*(voxel@[1,2,3])+3,atol=1e-5)
    with pytest.raises(ValueError,match='outside'):
        sample_sh(image,nib.affines.apply_affine(affine,[[-.1,0,0]]))


def test_vtk_world_roundtrip(tmp_path):
    points=np.array([[1.,2.,3.],[2.,2.,3.],[1.,3.,3.]])
    faces=np.array([[0,1,2]])
    path=tmp_path/'surface.vtk'
    write_surface(path,points,faces,{'roi':np.array([1,1,1],dtype=np.int32),'normals':np.tile([0.,0.,1.],(3,1))})
    p,f,d=read_surface(path)
    np.testing.assert_array_equal(p,points);np.testing.assert_array_equal(f,faces)
    np.testing.assert_array_equal(d['roi'],[1,1,1])
    header=path.read_bytes()
    assert header.startswith(b'# vtk DataFile Version 4.2\n')
    assert b'vtktypeint32' not in header
    assert b'roi 1 3 int\n' in header


@pytest.mark.parametrize('empty',[False,True])
def test_standard_tract_formats_and_empty_result(tmp_path,empty):
    affine=np.array([[1.25,0,0,-90],[0,1.25,0,-126],[0,0,1.25,-72],[0,0,0,1.]])
    image=nib.Nifti1Image(np.zeros((5,6,7)),affine)
    paths=[] if empty else [nib.affines.apply_affine(affine,np.array([[1.,1.,1.],[2.,3.,4.]]))]
    prefix=tmp_path/'tracks';write_tractograms(prefix,paths,image)
    assert validate_exports(prefix,paths,image)['streamlines']==len(paths)


def test_surface_validation_on_single_triangle(tmp_path):
    points=np.array([[1.,1.,1.],[2.,1.,1.],[1.,2.,1.]])
    geometry=build_geometry(points,np.array([[0,1,2]]))
    prepared=Prepared(geometry,np.ones((1,1,1),dtype=np.float32),np.ones(1,dtype=bool),np.zeros(1,dtype=np.int64),{})
    paths=[np.array([[1.2,1.2,1.],[1.5,1.3,1.]])]
    reference=nib.Nifti1Image(np.zeros((4,4,4)),np.eye(4))
    prefix=tmp_path/'single';write_tractograms(prefix,paths,reference)
    assert validate_exports(prefix,paths,reference,prepared)['surface_check']['max_plane_error_mm']<1e-7


def test_prepare_selects_one_gyrus_from_vtk_atlas(tmp_path):
    points=np.array([[1.,1.,1.],[2.,1.,1.],[1.,2.,1.],[5.,1.,1.],[6.,1.,1.],[5.,2.,1.]])
    faces=np.array([[0,1,2],[3,4,5]])
    surface=tmp_path/'atlas.vtk'
    write_surface(surface,points,faces,{'atlas':np.array([0,17,0,0,42,0])})
    image=tmp_path/'fod.nii'
    nib.save(nib.Nifti1Image(np.ones((8,4,4,1),dtype=np.float32),np.eye(4)),image)
    with pytest.raises(ValueError,match='multiple labels'):
        prepare(surface,image,shrink_mm=0,gyrus_array='atlas')
    p=prepare(surface,image,shrink_mm=0,gyrus_array='atlas',gyrus_label=17)
    np.testing.assert_array_equal(p.gyrus,[1,0])
    assert p.metadata['unified_gyrus_label']==17


def test_unified_gyrus_exports_validate_with_natural_stopping(tmp_path):
    points=np.array([[1.,1.,1.],[5.,1.,1.],[1.,5.,1.],[5.,5.,1.]])
    geometry=build_geometry(points,np.array([[0,1,3],[0,3,2]]))
    p=Prepared(geometry,np.ones((2,1,1),dtype=np.float32),np.ones(2,dtype=bool),
        np.zeros(2,dtype=np.int64),{},np.ones(2,dtype=np.int64))
    paths,report,d=track(p,TrackingConfig(n_seeds=20,threads=1,endpoint_policy='within_gyrus',
        stop_at_targets=False,max_turn_deg_per_mm=0.,min_length_mm=.1,max_length_mm=20.))
    assert len(paths)>0 and 'gyrus' not in report['stop_counts']
    reference=nib.Nifti1Image(np.zeros((7,7,3)),np.eye(4))
    prefix=tmp_path/'within';write_tractograms(prefix,paths,reference)
    validation=validate_exports(prefix,paths,reference,p,d)
    assert validation['all_endpoints_in_unified_gyrus']
