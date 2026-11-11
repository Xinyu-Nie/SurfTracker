import itertools
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse import csgraph

from swmtrack.open_preprocessing import (mesh_operators,curvature_feature,weighted_graph,
    potts_cut,vertex_links,thin_region,hamilton_jacobi_flux,thin_region_flux,
    skeleton_branches,prune_skeleton,gaussian_costs)
from swmtrack.curvature_estimators import (quadric_principal_curvatures,
                                           implicit_grid_principal_curvatures)


def sphere(refinements=3):
    p=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.],[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]])
    f=np.array([[0,1,2],[0,2,3],[0,3,4],[0,4,1],[5,2,1],[5,3,2],[5,4,3],[5,1,4]])
    for _ in range(refinements):
        points=list(p);midpoints={};faces=[]
        def midpoint(a,b):
            key=tuple(sorted((a,b)))
            if key not in midpoints:
                value=(p[a]+p[b])/2;value/=np.linalg.norm(value)
                midpoints[key]=len(points);points.append(value)
            return midpoints[key]
        for a,b,c in f:
            ab,bc,ca=midpoint(a,b),midpoint(b,c),midpoint(c,a)
            faces.extend([(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)])
        p=np.array(points);f=np.array(faces)
    return p,f


def test_cotan_fem_preserves_constants_area_and_rigid_frames():
    p,f=sphere(2);lap,mass,_=mesh_operators(p,f)
    np.testing.assert_allclose(lap@np.ones(len(p)),0,atol=1e-14)
    np.testing.assert_allclose(lap.toarray(),lap.toarray().T,atol=1e-14)
    area=np.linalg.norm(np.cross(p[f[:,1]]-p[f[:,0]],p[f[:,2]]-p[f[:,0]]),axis=1).sum()/2
    assert mass.sum()==pytest.approx(area)
    rotation=np.array([[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]])
    moved=p@rotation.T+[12.,-7.,3.]
    _,a=curvature_feature(p,f,0,0);_,b=curvature_feature(moved,f,0,0)
    np.testing.assert_allclose(a['mean_curvature'],b['mean_curvature'],atol=1e-12)
    np.testing.assert_allclose(a['gaussian_curvature'],b['gaussian_curvature'],atol=1e-12)


def test_curvature_has_known_sphere_sign_scaling_and_gauss_bonnet():
    p,f=sphere(3);_,values=curvature_feature(p,f,0,0)
    _,mass,_=mesh_operators(p,f)
    assert np.median(values['mean_curvature'])==pytest.approx(-1,rel=.03)
    assert np.median(values['gaussian_curvature'])==pytest.approx(1,rel=.03)
    assert np.sum(values['gaussian_curvature']*mass)==pytest.approx(4*np.pi,abs=1e-10)
    _,scaled=curvature_feature(3*p,f,0,0)
    np.testing.assert_allclose(scaled['mean_curvature'],values['mean_curvature']/3,atol=1e-12)
    np.testing.assert_allclose(scaled['gaussian_curvature'],values['gaussian_curvature']/9,atol=1e-12)


def test_cylinder_distinguishes_principal_order_and_concave_convex_folds():
    # A sphere cannot expose a swapped principal branch: its eigenvalues agree.
    # A unit cylinder must have curvatures (-1, 0) in the outward-negative
    # convention, and (0, 1) when viewed as a concave cylindrical wall.
    count=64;layers=9
    angle=np.arange(count)*2*np.pi/count
    p=np.array([[np.cos(t),np.sin(t),z] for z in np.linspace(-2,2,layers) for t in angle])
    f=[]
    for z in range(layers-1):
        for i in range(count):
            a=z*count+i;b=z*count+(i+1)%count;c=a+count;d=b+count
            f.extend([(a,b,c),(b,d,c)])
    f=np.array(f);interior=np.arange(count, count*(layers-1))
    convex,c=curvature_feature(p,f,0,0)
    concave,q=curvature_feature(p,f[:,[0,2,1]],0,0)
    np.testing.assert_allclose(convex[interior],-1,atol=.002)
    np.testing.assert_allclose(concave[interior],0,atol=1e-10)
    np.testing.assert_allclose(c['principal_curvature_min'][interior],-q['principal_curvature_max'][interior],atol=1e-10)


@pytest.mark.parametrize('method', ['quadric', 'implicit-grid'])
def test_alternative_estimators_recover_sphere_scale_and_sign(method):
    p, f = sphere(3)
    p *= 10
    if method == 'quadric':
        low, high, details = quadric_principal_curvatures(p, f, radius_mm=5)
    else:
        low, high, details = implicit_grid_principal_curvatures(p, f,
            spacing_mm=1, derivative_sigma_mm=1, band_mm=4)
        assert details['candidate_spotcheck_max_error_mm'] < .05
        assert details['gradient_magnitude_min'] > .8
    assert np.median(low) == pytest.approx(-.1, abs=.015)
    assert np.median(high) == pytest.approx(-.1, abs=.015)
    assert np.all(low <= high)
    assert np.isfinite(low).all() and np.isfinite(high).all()


def test_grid_curvature_uses_unit_tangent_normal_below_old_fallback(monkeypatch):
    import swmtrack.curvature_estimators as estimators
    p, f = sphere(3)
    p *= 10
    axis = np.arange(-14., 15.)
    x, y, z = np.meshgrid(axis, axis, axis, indexing='ij')
    # An attenuated signed-distance field has gradient magnitude near 0.3.
    # The tangent basis must stay orthonormal even though the curvature
    # denominator is floored at 0.5.
    field = .3 * (np.sqrt(x*x + y*y + z*z) - 10)
    monkeypatch.setattr(estimators, '_surface_signed_distance_grid',
                        lambda *args, **kwargs: (field.astype(np.float32),
                                                  np.array([-14., -14., -14.]), {}))
    low, high, details = estimators.implicit_grid_principal_curvatures(
        p, f, gradient_fallback_threshold=.2)
    assert details['quadric_fallback_vertices'] == 0
    assert details['gradient_magnitude_median'] == pytest.approx(.3, abs=.025)
    assert np.median(low) == pytest.approx(-.06, abs=.01)
    assert np.median(high) == pytest.approx(-.06, abs=.01)


@pytest.mark.parametrize('method', ['quadric', 'implicit-grid'])
def test_alternative_estimators_distinguish_cylinder_principal_branches(method):
    count, layers, radius = 64, 9, 10
    angles = np.arange(count) * 2*np.pi/count
    p = np.array([[radius*np.cos(t), radius*np.sin(t), z]
                  for z in np.linspace(-5, 5, layers) for t in angles])
    faces = []
    for layer in range(layers-1):
        for i in range(count):
            a=layer*count+i; b=layer*count+(i+1)%count; c=a+count; d=b+count
            faces.extend([(a,b,c),(b,d,c)])
    bottom,top=len(p),len(p)+1
    p=np.vstack([p,[0,0,-5],[0,0,5]])
    for i in range(count):
        faces.extend([(bottom,(i+1)%count,i),
                      (top,(layers-1)*count+i,(layers-1)*count+(i+1)%count)])
    faces=np.array(faces)
    middle=np.arange(3*count,6*count)
    if method == 'quadric':
        low,high,_=quadric_principal_curvatures(p,faces,radius_mm=3,source_ids=middle)
    else:
        low,high,_=implicit_grid_principal_curvatures(p,faces)
        low,high=low[middle],high[middle]
    assert np.median(low)==pytest.approx(-.1,abs=.015)
    assert np.median(high)==pytest.approx(0,abs=.02)


def test_quadric_expands_sparse_connected_neighborhoods():
    p, f = sphere(3)
    p *= 10
    low, high, details = quadric_principal_curvatures(p, f, radius_mm=3)
    assert details['expanded_neighborhood_count'] > 0
    assert details['maximum_effective_radius_mm'] > 3
    assert np.isfinite(low).all() and np.isfinite(high).all()


@pytest.mark.parametrize('method', ['quadric', 'implicit-grid'])
def test_alternative_feature_keeps_original_surface_and_vertex_ids(tmp_path, method):
    from swmtrack.preprocessing import preprocess, PreprocessConfig
    from swmtrack.io import write_surface, read_surface
    p, f = sphere(3)
    p *= 10
    source = tmp_path / 'input.vtk'
    write_surface(source, p, f, {'aparc': np.full(len(p), 17),
                                 'source_vertex_id': np.arange(len(p))})
    config = PreprocessConfig(curvature_method=method, mesh_smoothing=0,
                              quadric_radius_mm=5)
    meta = preprocess(source, tmp_path / 'output', config=config, stop_after='features')
    saved, saved_faces, fields = read_surface(tmp_path / 'output/preprocessed.world.vtk')
    np.testing.assert_array_equal(saved, p)
    np.testing.assert_array_equal(saved_faces, f)
    np.testing.assert_array_equal(fields['source_vertex_id'], np.arange(len(p)))
    assert np.isfinite(fields['k1']).all()
    assert meta['config']['curvature_method'] == method


def test_previous_method_segmentation_cache_is_rejected(tmp_path):
    from swmtrack.preprocessing import preprocess,PreprocessConfig
    from swmtrack.io import write_surface,sha256
    import json
    p,f=sphere(1);source=tmp_path/'input.vtk'
    write_surface(source,p,f,{'aparc':np.full(len(p),17)})
    cache=tmp_path/'cache';cache.mkdir()
    (cache/'preprocessing.json').write_text(json.dumps({'surface_sha256':sha256(source),
        'atlas_array':'aparc','method_version':'mesh-cotan-potts-link-thinning-v1'}))
    with pytest.raises(ValueError,match='method version differs'):
        preprocess(source,tmp_path/'out',reuse_segmentation=cache,stop_after='segmentation')


@pytest.mark.parametrize('seed',range(5))
def test_potts_cut_agrees_with_exhaustive_binary_energy(seed):
    rng=np.random.default_rng(seed);n=6
    costs=rng.integers(0,50,(n,2))/10
    edges=np.array([(i,i+1) for i in range(n-1)]+[(0,3),(2,5)])
    weights=rng.integers(0,15,len(edges))/10
    labels,report=potts_cut(costs,edges,weights)
    def energy(v):return costs[np.arange(n),v].sum()+weights[v[edges[:,0]]!=v[edges[:,1]]].sum()
    optimum=min(energy(np.array(v)) for v in itertools.product([0,1],repeat=n))
    assert energy(labels)==pytest.approx(optimum,abs=1e-10)
    assert report['floating_energy']==pytest.approx(optimum)


def test_graph_distance_uses_mesh_edges_and_thinning_preserves_disk_topology():
    p,f=sphere(3);graph,edges,lengths=weighted_graph(p,f)
    np.testing.assert_allclose(graph[edges[:,0],edges[:,1]].A.ravel(),lengths)
    mask=p[:,2]>-.1
    distance=csgraph.dijkstra(graph,indices=np.flatnonzero(~mask),min_only=True,directed=False)
    indptr,links=vertex_links(f,len(p));result=thin_region(mask,distance,indptr,links)
    assert 1<result.sum()<mask.sum()/2
    assert np.all(mask[result])
    for region in (mask,~mask,result,~result):
        ids=np.flatnonzero(region)
        assert csgraph.connected_components(graph[ids][:,ids],directed=False)[0]==1
    branches=skeleton_branches(graph,result)
    for b in branches:
        assert np.all(graph[b[:-1],b[1:]].A.ravel()>0)


def test_hj_flux_thinning_preserves_a_connected_surface_region():
    p,f=sphere(3)
    graph,_,_=weighted_graph(p,f)
    mask=p[:,2]>-.1
    distance,flux=hamilton_jacobi_flux(p,f,graph,mask)
    assert np.isfinite(distance).all() and np.isfinite(flux).all()
    assert np.quantile(flux[mask],.05)<np.median(flux[mask])
    indptr,links=vertex_links(f,len(p))
    result=thin_region_flux(mask,flux,float(np.percentile(flux[mask],5)),indptr,links)
    assert 0<result.sum()<mask.sum()/2
    assert np.all(mask[result])
    ids=np.flatnonzero(result)
    assert csgraph.connected_components(graph[ids][:,ids],directed=False)[0]==1


def test_area_guard_preserves_large_sulcus_when_skeleton_misses_it():
    from swmtrack.preprocessing import reconstruct_regions
    row=np.arange(5)
    graph=sparse.coo_matrix((np.ones(10),(np.r_[row,row+1],np.r_[row+1,row])),
                            shape=(6,6)).tocsr()
    initial=np.array([0,0,1,1,0,0])
    result=reconstruct_regions(initial,[],[np.array([0])],graph,
        valid=np.ones(6,dtype=bool),vertex_area=np.array([1,1,10,10,1,1]),
        max_unsupported_area_mm2=5)
    np.testing.assert_array_equal(result,[0,0,1,1,1,1])


def test_hj_feedback_runs_without_native_tools_or_moving_surface(tmp_path):
    from swmtrack.preprocessing import preprocess,PreprocessConfig
    from swmtrack.io import write_surface,read_surface
    p,f=sphere(3)
    p*=10
    source=tmp_path/'source.vtk';feature=tmp_path/'feature.raw'
    write_surface(source,p,f,{'aparc':np.full(len(p),17),
                              'source_vertex_id':np.arange(len(p))})
    (p[:,2]/10+.2*p[:,0]/10).astype('<f4').tofile(feature)
    meta=preprocess(source,tmp_path/'out',config=PreprocessConfig(
        skeleton_method='hj-flux',skeleton_relabeling='pruned-branches',
        gyral_pruning_mm=0,sulcal_pruning_mm=0),input_feature=feature,
        stop_after='segmentation')
    saved,saved_faces,arrays=read_surface(tmp_path/'out/preprocessed.world.vtk')
    np.testing.assert_array_equal(saved,p)
    np.testing.assert_array_equal(saved_faces,f)
    np.testing.assert_array_equal(arrays['source_vertex_id'],np.arange(len(p)))
    assert set(arrays['sulcal_segmentation'])=={0,1}
    assert np.isfinite(arrays['sulcal_flux']).all()
    assert meta['method_details']['skeleton_relabeling']['method']=='pruned-branches'


def test_physical_leaf_pruning_keeps_long_branch_and_removes_short_spur():
    p=np.array([[0.,0.,0.],[1.,0.,0.],[2.,0.,0.],[1.,.2,0.]])
    row=np.array([0,1,1,2,1,3]);col=np.array([1,0,2,1,3,1])
    graph=sparse.coo_matrix((np.linalg.norm(p[row]-p[col],axis=1),(row,col)),shape=(4,4)).tocsr()
    branches,_=prune_skeleton(p,graph,np.ones(4,dtype=bool),.5)
    assert len(branches)==1
    assert set(branches[0])=={0,1,2}


def test_cycles_are_traced_without_inventing_edges():
    graph=sparse.csr_matrix([[0.,1.,1.],[1.,0.,1.],[1.,1.,0.]])
    branches=skeleton_branches(graph,np.ones(3,dtype=bool))
    assert len(branches)==1 and len(branches[0])==4
    assert branches[0][0]==branches[0][-1]


def test_gaussian_prototypes_resist_tiny_extreme_outlier_class():
    rng=np.random.default_rng(62)
    signal=np.r_[rng.normal(-1,.2,500),rng.normal(1,.2,500),100000.]
    costs,fit=gaussian_costs(signal,np.ones(len(signal)),np.ones(len(signal),dtype=bool))
    assert np.isfinite(costs).all()
    assert np.mean(np.argmin(costs[:500],axis=1)==0)>.99
    assert np.mean(np.argmin(costs[500:1000],axis=1)==1)>.99
    assert -1.2<fit['means'][0]<-.8 and .8<fit['means'][1]<1.2
    for invalid in (np.ones(10),np.r_[np.ones(999),2.]):
        with pytest.raises(ValueError,match='Curvature'):
            gaussian_costs(invalid,np.ones(len(invalid)),np.ones(len(invalid),dtype=bool))


def test_preprocessing_runs_with_no_legacy_payload(tmp_path,monkeypatch):
    from swmtrack import native
    from swmtrack.preprocessing import preprocess,PreprocessConfig
    from swmtrack.io import write_surface,read_surface
    def forbidden(*args,**kwargs):raise AssertionError('Legacy executable resolution was attempted')
    monkeypatch.setattr(native,'bundle_root',forbidden)
    p,f=sphere(3)
    # Fold a closed surface to exercise real geometry, segmentation and thinning.
    p*=np.where(p[:,2]>.2,.7,1.)[:,None]
    source=tmp_path/'input.vtk'
    write_surface(source,p,f,{'aparc':np.full(len(p),17),'source_vertex_id':np.arange(len(p))})
    meta=preprocess(source,tmp_path/'out',stop_after='segmentation',
                    config=PreprocessConfig(gyral_pruning_mm=0,sulcal_pruning_mm=0))
    actual,faces,arrays=read_surface(tmp_path/'out/preprocessed.world.vtk')
    np.testing.assert_array_equal(actual,p);np.testing.assert_array_equal(faces,f)
    np.testing.assert_array_equal(arrays['source_vertex_id'],np.arange(len(p)))
    assert set(arrays['sulcal_segmentation'])=={0,1}
    assert meta['native_calls']==[] and meta['native_tools']=={} and meta['backend']=='python'
    write_surface(source,p,f[:,[0,2,1]],{'aparc':np.full(len(p),17)})
    with pytest.raises(ValueError,match='Outward'):
        preprocess(source,tmp_path/'inward',stop_after='features')
