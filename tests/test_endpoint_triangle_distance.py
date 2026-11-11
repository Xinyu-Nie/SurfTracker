import numpy as np
import pytest
from swmtrack.geometry import build_geometry
from swmtrack.validation import point_face_error, point_triangle_distance


def test_closed_triangle_distance_interior_edges_vertices_and_normal_offsets():
    g = build_geometry(np.array([[0.,0.,0.],[2.,0.,0.],[0.,2.,0.]]), np.array([[0,1,2]]))
    points = np.array([[.5,.5,.3],[1.5,1.5,0.],[-1.,-1.,0.],[.5,-.4,.3],[1.,1.,0.]])
    result = point_triangle_distance(g, points, np.zeros(len(points),int))
    np.testing.assert_allclose(result,[.3, np.sqrt(.5), np.sqrt(2), .5, 0.],atol=1e-14)


def test_float32_endpoint_on_thin_triangle_uses_physical_distance():
    # Regression from a valid on-edge endpoint in 101107 RH. These small
    # coordinate errors exceed the old relative barycentric threshold.
    tri = np.array([[40.72160496608995,-55.50548395205662,29.74297588536728],
                    [40.523606865764215,-55.220028692567965,28.684182937416153],
                    [40.63603765055644,-55.39226246092323,29.31205401374409]])
    g = build_geometry(tri,np.array([[0,1,2]]))
    point = np.array([[40.560516357421875,-55.27323913574219,28.88155174255371]],np.float32)
    _, bary = point_face_error(g,point,np.array([0]))
    assert bary[0] < -2e-4
    distance = point_triangle_distance(g,point,np.array([0]))[0]
    assert distance == pytest.approx(1.3293360272892758e-6,rel=1e-7)
    assert distance < 1e-4
    displaced = point.astype(float)+.001*g.frames[0,:,2]
    assert point_triangle_distance(g,displaced,np.array([0]))[0] > 9e-4


def test_distance_checks_selected_triangle_and_does_not_accept_the_whole_plane():
    points = np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],
                       [4.,0.,0.],[5.,0.,0.],[4.,1.,0.]])
    g = build_geometry(points,np.array([[0,1,2],[3,4,5]]))
    query = np.array([[4.2,.2,0.],[4.2,.2,0.]])
    result = point_triangle_distance(g,query,np.array([0,1]))
    assert result[0] > 3
    assert result[1] == 0
