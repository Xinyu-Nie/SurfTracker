import numpy as np
import pytest
from swmtrack.geometry import build_geometry
from swmtrack.tracking import Prepared, TrackingConfig, track, generate_seeds, unified_gyrus_faces, save_prepared, load_prepared
from swmtrack import ALGORITHM_VERSION
from swmtrack.kernel import categorical, random_uniform, ray_exit, sample_direction, trace_half, circular_value


def grid(n=4, size=12.):
    x, y = np.meshgrid(np.linspace(0, size, n+1), np.linspace(0, size, n+1))
    points = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
    faces = []
    for j in range(n):
        for i in range(n):
            a = j*(n+1)+i; b = a+1; c = a+n+1; d = c+1
            faces.extend([(a, b, d), (a, d, c)])
    return build_geometry(points, np.array(faces))


def prepared_grid(n=4):
    g = grid(n)
    coefficients = np.ones((len(g.faces), 1, 1), dtype=np.float32)
    targets = np.zeros(len(g.faces), dtype=np.int64)
    centers = g.points[g.faces].mean(axis=1)
    targets[centers[:, 0] < 2] = 1; targets[centers[:, 0] > 10] = 2
    seeds = (centers[:, 0] > 4) & (centers[:, 0] < 8)
    return Prepared(g, coefficients, seeds, targets, {})


def test_categorical_distribution_and_empty_support():
    weights = np.array([.1, .2, .7])
    counts = np.zeros(3)
    state = np.uint64(721)
    for _ in range(50000):
        k, state = categorical(weights, np.uint64(state)); counts[k] += 1
    np.testing.assert_allclose(counts/counts.sum(), weights, atol=.008)
    for w in (np.zeros(3), np.array([np.nan, -1., 0.])):
        assert categorical(w, np.uint64(1))[0] == -1


def test_flat_transport_and_ray_exit():
    g = grid()
    for f, e in zip(*np.where(g.neighbors >= 0)):
        matrix = g.transport[f, e]
        np.testing.assert_allclose(g.frames[g.neighbors[f, e], :, :2]@matrix,
                                   g.frames[f, :, :2], atol=1e-12)
    distance, edge, vertex = ray_exit(np.array([[0., 0.], [2., 0.], [0., 2.]]), .3, .3, 1., 0.)
    assert distance == pytest.approx(1.4)
    assert edge == 1 and not vertex


def test_folded_transport_preserves_shared_edge_and_tangent_norm():
    g = build_geometry(np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]), np.array([[0,1,2], [1,0,3]]))
    f, e = 0, 0
    matrix = g.transport[f,e]
    np.testing.assert_allclose(matrix.T@matrix, np.eye(2), atol=1e-12)
    edge = np.array([1.,0.,0.])
    np.testing.assert_allclose(g.frames[1,:,:2]@matrix@(g.frames[0,:,:2].T@edge), edge, atol=1e-12)


def test_invalid_meshes_raise():
    with pytest.raises(ValueError, match="Degenerate"):
        build_geometry(np.array([[0.,0.,0.], [1.,0.,0.], [2.,0.,0.]]), np.array([[0,1,2]]))
    with pytest.raises(ValueError, match="winding"):
        build_geometry(np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]), np.array([[0,1,2],[0,1,3]]))


def test_zero_or_below_threshold_fod_stops_without_fallback():
    g = grid(); f = 0; xy = g.xy[f]; p = xy.mean(axis=0)
    for value in (0., .001):
        result = sample_direction(np.array([[value]]), xy, p[0], p[1], 0., np.pi,
                                  np.pi/90, .01, np.uint64(9))
        assert not result[2]


def test_fixed_mm_straight_path_survives_mesh_refinement():
    end_points, lengths = [], []
    for n in (2, 8):
        g = grid(n)
        seed = np.array([2.31, 3.47, 0.])
        f = -1
        for k, triangle in enumerate(g.points[g.faces]):
            bary = np.linalg.lstsq(np.vstack((triangle.T, np.ones(3))), np.r_[seed,1], rcond=None)[0]
            if np.min(bary) > 1e-8:
                f = k; break
        assert f >= 0
        world_direction = np.array([1., .17, 0.]); world_direction /= np.linalg.norm(world_direction)
        d = g.frames[f,:,:2].T@world_direction
        coeff = np.ones((len(g.faces),1,1), dtype=np.float32)
        path,count,status,_,length,*_ = trace_half(g.points,g.faces,g.frames,g.xy,g.neighbors,g.transport,
            coeff,np.zeros(len(g.faces),dtype=np.int64),f,seed,d[0],d[1],.37,0.,np.pi/90,.01,30.,1000,False,np.uint64(4))
        assert status == 3
        end_points.append(path[count-1]); lengths.append(length)
    np.testing.assert_allclose(end_points[0], end_points[1], atol=2e-10)
    np.testing.assert_allclose(lengths[0], lengths[1], atol=2e-10)


def test_thread_and_batch_reproducibility_and_endpoint_policy():
    p = prepared_grid(8)
    common = dict(n_seeds=60,seed_mode="area",min_length_mm=0.1,max_length_mm=40.,step_mm=.5)
    a,ra,da = track(p, TrackingConfig(**common,threads=1,batch_size=17))
    b,rb,db = track(p, TrackingConfig(**common,threads=2,batch_size=23))
    assert ra['accepted'] > 0
    assert len(a) == len(b)
    for x,y in zip(a,b): np.testing.assert_array_equal(x,y)
    np.testing.assert_array_equal(da['stops'],db['stops'])
    accepted = da['stops'][da['accepted']]
    assert np.all(np.sort(accepted,axis=1)==[1,2])
    _,ru,_ = track(p, TrackingConfig(**common,threads=2,endpoint_policy="unrestricted"))
    assert ru['stop_at_targets'] is False
    assert not ('gyrus_a' in ru['stop_counts'] or 'gyrus_b' in ru['stop_counts'])


def test_empty_tracks_are_valid():
    p=prepared_grid();p.coefficients[:]=0
    paths,report,d=track(p,TrackingConfig(n_seeds=4,threads=1,endpoint_policy="unrestricted"))
    assert paths==[] and report['accepted']==0
    assert report['stop_counts']=={'no_supported_direction':8}


def test_worker_pool_limit_preserves_seeded_paths(monkeypatch):
    import numba
    p = prepared_grid(8)
    common = dict(n_seeds=40, min_length_mm=.1, max_length_mm=40.)
    reference, _, _ = track(p, TrackingConfig(**common, threads=1))
    original_threads = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        with monkeypatch.context() as patch:
            patch.setattr(numba.config, 'NUMBA_NUM_THREADS', 1)
            paths, report, _ = track(p, TrackingConfig(**common, threads=4))
    finally:
        numba.set_num_threads(original_threads)
    assert report['config']['threads'] == 4
    assert report['effective_threads'] == report['thread_pool_limit'] == 1
    assert len(paths) == len(reference) > 0
    for actual, expected in zip(paths, reference):
        np.testing.assert_array_equal(actual, expected)


def test_area_seeds_follow_physical_area_and_uniform_barycentric_positions():
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[5.,0.,0.],[7.,0.,0.],[5.,2.,0.]])
    g=build_geometry(points,np.array([[0,1,2],[3,4,5]]))
    p=Prepared(g,np.ones((2,1,1),dtype=np.float32),np.ones(2,dtype=bool),np.zeros(2,dtype=np.int64),{})
    faces,seeds=generate_seeds(p,TrackingConfig(n_seeds=30000,seed_mode='area'))
    assert np.mean(faces==0)==pytest.approx(.2,abs=.008)
    for face in (0,1):
        np.testing.assert_allclose(seeds[faces==face].mean(axis=0),g.points[g.faces[face]].mean(axis=0),atol=.025)


def test_exact_folded_path_matches_unfolded_geometry():
    points=np.array([[0.,0.,0.],[5.,0.,0.],[5.,0.,5.],[0.,10.,0.],[5.,10.,0.],[5.,10.,5.]])
    faces=np.array([[0,1,4],[0,4,3],[1,2,5],[1,5,4]])
    g=build_geometry(points,faces)
    seed=np.array([1.3,2.7,0.]);direction=np.array([1.,.1,0.]);direction/=np.linalg.norm(direction)
    f=1;d=g.frames[f,:,:2].T@direction
    result=trace_half(g.points,g.faces,g.frames,g.xy,g.neighbors,g.transport,
        np.ones((4,1,1),dtype=np.float32),np.zeros(4,dtype=np.int64),f,seed,d[0],d[1],
        .37,0.,np.pi/90,.01,50.,1000,False,np.uint64(12))
    path,count,status,_,length,*_=result
    assert status==3
    np.testing.assert_allclose(path[count-1],[5.,2.7+.1*(10-1.3),5.],atol=1e-10)
    assert length==pytest.approx((10-1.3)*np.sqrt(1.01),abs=1e-10)


def test_overlapping_targets_do_not_take_gyrus_precedence():
    p=prepared_grid();p.targets[:]=3
    paths,report,d=track(p,TrackingConfig(n_seeds=10,threads=1,min_length_mm=.01,max_length_mm=40.))
    assert paths==[]
    assert report['stop_counts'].get('overlapping_targets',0)>0
    assert not ('gyrus_a' in report['stop_counts'] or 'gyrus_b' in report['stop_counts'])


def test_post_projection_clipping_is_applied_at_evaluation():
    # g(theta)=0.1+0.2*cos(2theta), including a genuinely negative sector.
    coefficients=np.array([.1,.2,0.])
    assert circular_value(coefficients,0.)==pytest.approx(.3)
    assert circular_value(coefficients,np.pi/2)==0.


def test_nonmanifold_edge_is_rejected():
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,-1.,0.]])
    with pytest.raises(ValueError,match='Nonmanifold'):
        build_geometry(points,np.array([[0,1,2],[1,0,3],[0,1,4]]))


def test_within_gyrus_requires_both_endpoints_in_one_mask():
    p=prepared_grid(8)
    # Both marked ends belong to ONE synthetic gyrus G; A/B masks are absent.
    p.gyrus=(p.targets>0).astype(np.int64)
    p.targets[:]=0
    paths,report,d=track(p,TrackingConfig(n_seeds=120,threads=2,min_length_mm=.1,
        max_length_mm=40.,endpoint_policy='within_gyrus'))
    assert len(paths)>0 and report['stop_at_targets']
    assert report['stop_counts'].get('gyrus',0)>0
    assert 'gyrus_a' not in report['stop_counts'] and 'gyrus_b' not in report['stop_counts']
    ids=d['accepted_trial_ids']
    assert np.all(d['stops'][ids]==1)
    assert np.all(p.gyrus[d['endpoint_faces'][ids]]==1)


def test_within_gyrus_rejects_only_one_endpoint_touching():
    p=prepared_grid(8)
    p.gyrus=(p.targets==1).astype(np.int64)
    # Straight lines from the center can reach the left-hand mask at only one end.
    paths,report,d=track(p,TrackingConfig(n_seeds=100,threads=1,min_length_mm=.1,
        max_length_mm=40.,max_turn_deg_per_mm=0.,endpoint_policy='within_gyrus'))
    assert paths==[]
    assert np.any(d['stops']==1) and np.any(d['stops']==3)


def test_within_gyrus_enforces_endpoints_without_stopping_on_target_entry():
    assert TrackingConfig(endpoint_policy='within_gyrus').validate() is True
    assert TrackingConfig(endpoint_policy='within_gyrus',stop_at_targets=False).validate() is False
    p=prepared_grid(8);p.gyrus=np.ones(len(p.geometry.faces),dtype=np.int64)
    paths,report,d=track(p,TrackingConfig(n_seeds=40,threads=2,min_length_mm=.1,
        max_length_mm=40.,max_turn_deg_per_mm=0.,endpoint_policy='within_gyrus',stop_at_targets=False))
    assert len(paths)>0 and not report['stop_at_targets']
    assert not np.any(d['stops']==1)
    assert np.all(d['endpoint_labels'][d['accepted']]==1)


def test_within_gyrus_requires_an_explicit_nonempty_single_label():
    with pytest.raises(ValueError,match='one nonempty gyrus'):
        track(prepared_grid(),TrackingConfig(endpoint_policy='within_gyrus'))


def test_between_gyri_enforces_endpoints_without_stopping_on_target_entry():
    p=prepared_grid(8)
    paths,report,d=track(p,TrackingConfig(n_seeds=100,threads=2,min_length_mm=.1,
        max_length_mm=40.,max_turn_deg_per_mm=0.,endpoint_policy='between_gyri',stop_at_targets=False))
    assert len(paths)>0 and not report['stop_at_targets']
    assert not np.any(np.isin(d['stops'],[1,2]))
    assert np.all(np.sort(d['endpoint_labels'][d['accepted']],axis=1)==[1,2])


def test_unified_gyrus_does_not_silently_merge_atlas_labels():
    faces=np.array([[0,1,2],[3,4,5]])
    labels=np.array([0,17,0,0,42,0])
    with pytest.raises(ValueError,match='multiple labels'):
        unified_gyrus_faces(labels,faces)
    np.testing.assert_array_equal(unified_gyrus_faces(labels,faces,label=17),[1,0])
    np.testing.assert_array_equal(unified_gyrus_faces(labels,faces,label=42),[0,1])
    with pytest.raises(ValueError,match='positive'):
        unified_gyrus_faces(labels,faces,label=0)


def test_unified_gyrus_cache_roundtrip_and_old_cache_compatibility(tmp_path):
    p=prepared_grid();p.metadata={'algorithm_version':ALGORITHM_VERSION}
    p.gyrus=(p.targets==1).astype(np.int64)
    path=tmp_path/'with_gyrus.npz';save_prepared(path,p)
    np.testing.assert_array_equal(load_prepared(path).gyrus,p.gyrus)
    p.gyrus=None
    path=tmp_path/'without_gyrus.npz';save_prepared(path,p)
    assert load_prepared(path).gyrus is None
    with pytest.raises(ValueError,match='one nonempty gyrus'):
        track(load_prepared(path),TrackingConfig(endpoint_policy='within_gyrus'))
