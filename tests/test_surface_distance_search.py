import numpy as np
from scipy.spatial import cKDTree
from swmtrack.geometry import build_geometry
from swmtrack.validation import surface_point_distances


def geometry_with_distant_containing_triangle_center():
    triangles=[np.array([[0.,0.,0.],[100.,0.,0.],[0.,100.,0.]])]
    for i in range(40):
        x=.01+i*.003
        triangles.append(np.array([[x,.01,.01],[x+.001,.01,.01],[x,.011,.01]]))
    points=np.vstack(triangles)
    return build_geometry(points,np.arange(len(points)).reshape(-1,3))


def test_bounded_fallback_finds_triangle_missed_by_nearest_32_centers():
    g=geometry_with_distant_containing_triangle_center();query=np.array([[.01,.01,0.]])
    tree=cKDTree(g.points[g.faces].mean(axis=1));_,ids=tree.query(query,k=32)
    assert 0 not in ids
    distances,faces,fallback=surface_point_distances(g,query)
    assert distances[0]==0 and faces[0]==0 and fallback==1


def test_bounded_fallback_does_not_accept_an_off_surface_point():
    g=geometry_with_distant_containing_triangle_center()
    distances,_,fallback=surface_point_distances(g,np.array([[.01,.01,.1]]))
    assert distances[0]>.08 and fallback==1


def test_single_triangle_query_preserves_plane_error_and_edge_distance():
    g=build_geometry(np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]),np.array([[0,1,2]]))
    distances,faces,_=surface_point_distances(g,np.array([[.2,.2,1e-6],[.8,.8,0.]]))
    np.testing.assert_allclose(distances,[1e-6,.6/np.sqrt(2)],atol=1e-14)
    np.testing.assert_array_equal(faces,[0,0])
