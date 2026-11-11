"""Independent mesh preprocessing using NumPy, SciPy and Numba.

Implements standard cotangent FEM/angle-defect formulas, a two-Gaussian Potts
model solved by SciPy's integer max-flow, edge-graph distances, and simple-point
thinning on triangular vertex links. No GW, Orion3D or legacy max-flow code is
used. This is a different estimator from the legacy implicit-curvature method.
"""

import heapq
import json
from pathlib import Path
import numpy as np
from numba import njit
from scipy import sparse
from scipy.sparse import csgraph
from scipy.sparse.linalg import splu

from .geometry import vertex_normals


METHOD_VERSION = 'python-consistent-mass-agreement-candidate-v6'


def check_python():
    """Check required distributions and the integer max-flow runtime."""
    from importlib.metadata import version
    versions={name:version(name) for name in ('numpy','scipy','numba','nibabel','meshio')}
    result=csgraph.maximum_flow(sparse.csr_matrix([[0,3],[0,0]],dtype=np.int64),0,1,method='dinic')
    if result.flow_value!=3:raise RuntimeError('SciPy maximum-flow runtime check failed')
    return {'backend':'python','method_version':METHOD_VERSION,'dependencies':versions,
            'legacy_executables_required':False,'passed':True}


def mesh_operators(points, faces):
    """Positive cotangent stiffness and barycentric lumped vertex mass."""
    tri=np.asarray(points)[faces]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    twice_area=np.linalg.norm(cross,axis=1)
    if np.any(twice_area<=1e-12):raise ValueError('Degenerate geometry in curvature estimation')
    angles=[];cot=[]
    for i in range(3):
        a=tri[:,(i+1)%3]-tri[:,i];b=tri[:,(i+2)%3]-tri[:,i]
        dot=np.einsum('ij,ij->i',a,b)
        angles.append(np.arctan2(twice_area,dot));cot.append(dot/twice_area)
    row=np.concatenate([faces[:,(i+1)%3] for i in range(3)])
    col=np.concatenate([faces[:,(i+2)%3] for i in range(3)])
    weights=np.concatenate(cot)/2
    edge=sparse.coo_matrix((np.r_[weights,weights],(np.r_[row,col],np.r_[col,row])),
                          shape=(len(points),len(points))).tocsr()
    stiffness=sparse.diags(np.asarray(edge.sum(axis=1)).ravel())-edge
    mass=np.bincount(faces.ravel(),weights=np.repeat(twice_area/6,3),minlength=len(points))
    angle_sum=np.bincount(faces.ravel(),weights=np.stack(angles,axis=1).ravel(),minlength=len(points))
    if np.any(mass<=0):raise ValueError('Nonpositive vertex area')
    return stiffness.tocsc(),mass,angle_sum


def smoothing_mass_matrix(points, faces, area, kind):
    """Linear-triangle FEM mass or its barycentric lumped diagonal.

    The local consistent matrix is A/12 * [[2,1,1],[1,2,1],[1,1,2]].
    It is assembled independently from triangle areas; no native matrix is used.
    """
    if kind == 'lumped':
        return sparse.diags(area).tocsc()
    if kind != 'consistent':
        raise ValueError('Smoothing mass must be lumped or consistent')
    tri=np.asarray(points)[faces]
    triangle_area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    if not np.isfinite(triangle_area).all() or np.any(triangle_area<=0):
        raise ValueError('Consistent mass requires finite, nondegenerate triangles')
    row=np.repeat(faces,3,axis=1).ravel();col=np.tile(faces,(1,3)).ravel()
    local=np.ones((3,3))+np.eye(3)
    return sparse.coo_matrix(((triangle_area[:,None,None]*local/12).ravel(),(row,col)),
                             shape=(len(points),len(points))).tocsc()


def curvature_feature(points, faces, mesh_smoothing, feature_smoothing, *, smoothing_mass='lumped'):
    stiffness,mass,_=mesh_operators(points,faces)
    matrix=smoothing_mass_matrix(points,faces,mass,smoothing_mass)
    smoothed=splu(matrix+mesh_smoothing*stiffness).solve(matrix@points)
    lap,area,angles=mesh_operators(smoothed,faces)
    normal=vertex_normals(smoothed,faces)
    # Outward sphere has negative H; concave sulcal folds have positive H.
    mean=-.5*np.einsum('ij,ij->i',lap@smoothed,normal)/area
    gaussian=(2*np.pi-angles)/area
    # Reversing the curvature sign swaps the principal-curvature ordering:
    # min(-k1, -k2) = -max(k1, k2). The legacy feature is the largest
    # convex-positive principal curvature. Its concave-positive counterpart
    # is therefore H - sqrt(H**2-K), NOT H + sqrt(H**2-K).
    delta=np.sqrt(np.maximum(mean**2-gaussian,0))
    raw=mean-delta
    feature=splu(matrix+feature_smoothing*stiffness).solve(matrix@raw)
    if not np.isfinite(feature).all():raise ValueError('Nonfinite curvature feature')
    return feature.astype(np.float32),dict(smoothed=smoothed,mean_curvature=mean,
        gaussian_curvature=gaussian,principal_curvature_min=raw,
        principal_curvature_max=mean+delta,k1_unfiltered=raw,k1=feature)


def curvature_feature_options(points, faces, config):
    """Use an explicit alternative curvature estimator on a temporary mesh.

    Principal curvatures come from one fitted shape operator. The final feature
    is still smoothed on the original indexed mesh, then the original mesh is
    used unchanged for segmentation, patches, and tracking.
    """
    from .curvature_estimators import (quadric_principal_curvatures,
                                       implicit_grid_principal_curvatures)
    stiffness, mass, _ = mesh_operators(points, faces)
    matrix=smoothing_mass_matrix(points,faces,mass,config.smoothing_mass)
    smoothed = splu(matrix+config.mesh_smoothing*stiffness).solve(matrix@points)
    if config.curvature_method == 'quadric':
        low, high, details = quadric_principal_curvatures(
            smoothed, faces, radius_mm=config.quadric_radius_mm)
    elif config.curvature_method == 'implicit-grid':
        low, high, details = implicit_grid_principal_curvatures(
            smoothed, faces, spacing_mm=config.implicit_grid_spacing_mm,
            derivative_sigma_mm=config.implicit_derivative_sigma_mm,
            gradient_fallback_threshold=config.implicit_gradient_fallback_threshold,
            band_mm=config.implicit_band_mm,
            candidate_faces=config.implicit_candidate_faces,
            max_grid_voxels=config.implicit_max_grid_voxels)
    else:
        raise ValueError(f'Unsupported curvature_method: {config.curvature_method}')
    feature=splu(matrix+config.feature_smoothing*stiffness).solve(matrix@low)
    if not np.isfinite(feature).all():
        raise ValueError('Nonfinite curvature feature')
    arrays=dict(smoothed=smoothed,mean_curvature=.5*(low+high),
                gaussian_curvature=low*high,principal_curvature_min=low,
                principal_curvature_max=high,k1_unfiltered=low,k1=feature.astype(np.float32))
    return feature.astype(np.float32),arrays,{'method':config.curvature_method,
        'geometry_smoothing_mass':config.smoothing_mass,'feature_smoothing_mass':config.smoothing_mass,
        'coordinate_rounding':'none',**details}


def weighted_graph(points,faces):
    edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    edges=np.unique(edges,axis=0)
    lengths=np.linalg.norm(points[edges[:,0]]-points[edges[:,1]],axis=1)
    if np.any(lengths<=0):raise ValueError('Nonpositive graph edge length')
    graph=sparse.coo_matrix((np.r_[lengths,lengths],
        (np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(points),len(points))).tocsr()
    graph.sort_indices()
    return graph,edges,lengths


def gaussian_costs(feature,area,valid):
    """Robust two-class Gaussian prototypes with equal prior probabilities.

    Fit lower/upper median groups after winsorizing the fitting signal at its
    1st/99th percentiles. Equal priors prevent a tiny high-curvature outlier
    component from swallowing the intended two anatomical folding classes.
    Component 1 is the higher-score sulcal class under our sign convention.
    """
    x=np.asarray(feature,dtype=float);values=x[valid];weights=area[valid]
    if len(values)<2 or not np.isfinite(values).all():
        raise ValueError('At least two finite cortical curvature values are required')
    if np.ptp(values)<1e-7:raise ValueError('Curvature is constant; sulcal/gyral classes are unsupported')
    bounds=np.quantile(values,[.01,.99]);clipped=np.clip(values,*bounds)
    center=float(np.average(clipped,weights=weights));scale=max(float(np.std(clipped)),1e-8)
    z=(clipped-center)/scale;median=float(np.median(z));groups=(z<=median,z>median)
    if any(not g.any() for g in groups):
        raise ValueError('Curvature has no robust two-class split after outlier clipping')
    means=np.array([np.average(z[g],weights=weights[g]) for g in groups])
    variance=np.array([max(float(np.average((z[g]-means[i])**2,weights=weights[g])),.02) for i,g in enumerate(groups)])
    mix=np.full(2,.5)
    full=(np.clip(x,*bounds)-center)/scale
    costs=-np.log(mix)[None,:]+.5*np.log(2*np.pi*variance)[None,:]+.5*(full[:,None]-means)**2/variance
    costs-=costs.min(axis=1)[:,None]
    costs=np.minimum(costs,20)*area[:,None]
    costs[~valid]=0
    return costs,{'weights':mix.tolist(),'means':(means*scale+center).tolist(),
                  'variances':(variance*scale**2).tolist(),'winsorization_bounds':bounds.tolist(),
                  'fit':'area-weighted lower/upper median Gaussian prototypes; equal priors'}


def robust_unit_costs(feature, area, valid):
    """Unweighted whole-mesh costs with the unchanged robust cortical fit.

    Preserve cortical floating-point costs from gaussian_costs/area exactly.
    Evaluate that fitted, winsorized model at unknown vertices as well; they
    remain excluded from fitting and retain the existing post-cut masking.
    """
    original,details=gaussian_costs(feature,area,valid)
    p=np.column_stack((details['weights'],details['means'],details['variances']))
    x=np.clip(feature,*details['winsorization_bounds'])
    costs=.5*np.log(2*np.pi*p[:,2])-np.log(p[:,0])+.5*(x[:,None]-p[:,1])**2/p[:,2]
    costs-=costs.min(axis=1)[:,None]
    costs=np.minimum(costs,20)
    costs[valid]=original[valid]/area[valid,None]
    return costs,{**details,'unary_weighting':'unweighted','evaluation_domain':'whole mesh',
                  'fit_domain':'cortex only','likelihood_cap':20}


def potts_cut(costs, edges, pairwise):
    """Exact cut for quantized energy; report the floating energy error bound.

    Keep the sum of all directed capacities below 2**30 to bound int32 flow and
    residual arithmetic across supported SciPy versions. Zero terminal offsets
    are removed before this function; nonnegative floating costs are required.
    """
    costs=np.asarray(costs,dtype=float);pairwise=np.asarray(pairwise,dtype=float)
    if costs.ndim!=2 or costs.shape[1]!=2 or pairwise.shape!=(len(edges),):raise ValueError('Invalid Potts arrays')
    if not np.isfinite(costs).all() or not np.isfinite(pairwise).all() or np.any(costs<0) or np.any(pairwise<0):
        raise ValueError('Potts costs must be finite and nonnegative')
    n=len(costs);total=float(costs.sum()+2*pairwise.sum())
    scale=min(100000.,(2**30-2*(costs.size+2*len(edges)))/max(total,1e-12))
    if scale<=0:raise ValueError('Graph exceeds integer-capacity safety bound')
    unary=np.rint(costs*scale).astype(np.int64);pairs=np.rint(pairwise*scale).astype(np.int64)
    row=np.r_[np.full(n,n),np.arange(n),edges[:,0],edges[:,1]]
    col=np.r_[np.arange(n),np.full(n,n+1),edges[:,1],edges[:,0]]
    capacities=np.r_[unary[:,0],unary[:,1],pairs,pairs]
    graph=sparse.coo_matrix((capacities,(row,col)),shape=(n+2,n+2)).tocsr();graph.eliminate_zeros()
    flow=csgraph.maximum_flow(graph,n,n+1,method='dinic')
    residual=(graph-flow.flow).tocsr();residual.data=(residual.data>0).astype(np.int64);residual.eliminate_zeros()
    reachable=csgraph.breadth_first_order(residual,n,directed=True,return_predecessors=False)
    labels=np.zeros(n+2,dtype=np.int32);labels[reachable]=1;labels=labels[:n]
    energy=float(costs[np.arange(n),labels].sum()+pairwise[labels[edges[:,0]]!=labels[edges[:,1]]].sum())
    return labels,{'solver':'scipy.sparse.csgraph.maximum_flow/dinic','capacity_scale':scale,
        'quantized_flow':int(flow.flow_value),'floating_energy':energy,
        'worst_case_energy_rounding_error':float((n+len(edges))*.5/scale),
        'zeroed_positive_unary_capacities':int(np.sum((costs>0)&(unary==0))),
        'zeroed_positive_pairwise_capacities':int(np.sum((pairwise>0)&(pairs==0)))}


@njit(cache=True)
def _order_links(indptr,left,right):
    ordered=np.empty_like(left)
    for v in range(len(indptr)-1):
        start,end=indptr[v],indptr[v+1]
        current=left[start]
        for k in range(start,end):
            ordered[k]=current
            for j in range(start,end):
                if left[j]==current:
                    current=right[j];break
            else:raise ValueError('Incomplete vertex link')
            if current==ordered[start] and k<end-1:raise ValueError('Disconnected vertex link')
        if current!=ordered[start]:raise ValueError('Open vertex link')
    return ordered


def vertex_links(faces,n):
    center=faces.ravel();left=np.roll(faces,-1,axis=1).ravel();right=np.roll(faces,-2,axis=1).ravel()
    order=np.lexsort((left,center));indptr=np.r_[0,np.cumsum(np.bincount(center,minlength=n))]
    return indptr,_order_links(indptr,left[order],right[order])


@njit(cache=True)
def thin_region(mask,distances,indptr,neighbors):
    """Distance-priority simple-point deletion, preserving graph endpoints.

    For a manifold triangular mesh a simple boundary vertex has exactly one
    foreground arc in its cyclic link. Removing it preserves region topology.
    """
    active=mask.copy();queued=np.zeros(len(mask),dtype=np.bool_)
    heap=[(0.,0)];heap.pop()
    for v in range(len(mask)):
        if active[v]:
            for j in range(indptr[v],indptr[v+1]):
                if not active[neighbors[j]]:
                    heapq.heappush(heap,(distances[v],v));queued[v]=True;break
    while heap:
        _,v=heapq.heappop(heap);queued[v]=False
        if not active[v]:continue
        start,end=indptr[v],indptr[v+1];count=0;transitions=0
        previous=active[neighbors[end-1]]
        for j in range(start,end):
            value=active[neighbors[j]];count+=int(value);transitions+=int(value!=previous);previous=value
        if count<=1 or transitions!=2:continue
        active[v]=False
        for j in range(start,end):
            u=neighbors[j]
            if active[u] and not queued[u]:
                heapq.heappush(heap,(distances[u],u));queued[u]=True
    return active


def hamilton_jacobi_flux(points,faces,graph,mask):
    """Approximate the paper's one-ring average flux from a surface distance field.

    Mesh-edge shortest paths approximate triangular fast marching here. Face
    gradients of the distance transform are averaged at vertices before the
    one-ring flux in equation (5) of Shi et al. (2008) is evaluated.
    """
    mask=np.asarray(mask,dtype=bool)
    boundary=mask & (np.asarray(graph.astype(bool)@(~mask).astype(np.int32)).ravel()>0)
    if not boundary.any():
        raise ValueError('HJ flux needs a nonempty object boundary')
    distance=csgraph.dijkstra(graph,indices=np.flatnonzero(boundary),
                             min_only=True,directed=False)
    distance[~mask]=0
    tri=np.asarray(points)[faces]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    double_area=np.linalg.norm(cross,axis=1)
    if np.any(double_area<=1e-12):
        raise ValueError('Degenerate triangle in HJ flux calculation')
    normal=cross/double_area[:,None]
    values=distance[faces]
    face_gradient=(values[:,0,None]*np.cross(normal,tri[:,2]-tri[:,1])+
                   values[:,1,None]*np.cross(normal,tri[:,0]-tri[:,2])+
                   values[:,2,None]*np.cross(normal,tri[:,1]-tri[:,0]))/double_area[:,None]
    gradient=np.zeros_like(points,dtype=np.float64)
    area=np.zeros(len(points),dtype=np.float64)
    for corner in range(3):
        np.add.at(gradient,faces[:,corner],face_gradient*double_area[:,None])
        np.add.at(area,faces[:,corner],double_area)
    gradient/=np.maximum(area[:,None],1e-12)
    norm=np.linalg.norm(gradient,axis=1)
    gradient[norm>1e-8]/=norm[norm>1e-8,None]
    row=np.repeat(np.arange(len(points)),np.diff(graph.indptr))
    col=graph.indices
    direction=points[col]-points[row]
    edge_length=np.linalg.norm(direction,axis=1)
    term=np.einsum('ij,ij->i',direction,gradient[col])/edge_length
    flux=np.bincount(row,weights=term,minlength=len(points))/np.maximum(np.diff(graph.indptr),1)
    flux[~mask]=0
    return distance,flux


@njit(cache=True)
def thin_region_flux(mask,flux,threshold,indptr,neighbors):
    """Flux-ordered simple-point thinning with an endpoint flux threshold."""
    active=mask.copy();queued=np.zeros(len(mask),dtype=np.bool_)
    heap=[(0.,0)];heap.pop()
    for v in range(len(mask)):
        if active[v]:
            for j in range(indptr[v],indptr[v+1]):
                if not active[neighbors[j]]:
                    heapq.heappush(heap,(-flux[v],v));queued[v]=True;break
    while heap:
        _,v=heapq.heappop(heap);queued[v]=False
        if not active[v]:continue
        start,end=indptr[v],indptr[v+1];count=0;transitions=0
        previous=active[neighbors[end-1]]
        for j in range(start,end):
            value=active[neighbors[j]];count+=int(value);transitions+=int(value!=previous);previous=value
        if count==0 or transitions!=2:continue
        if count==1 and flux[v]<=threshold:continue
        active[v]=False
        for j in range(start,end):
            u=neighbors[j]
            if active[u] and not queued[u]:
                heapq.heappush(heap,(-flux[u],u));queued[u]=True
    return active


def skeleton_branches(graph,mask):
    """Trace each undirected skeleton edge once, including closed cycles."""
    ids=np.flatnonzero(mask);sub=graph[ids][:,ids].tocsr();degree=np.diff(sub.indptr)
    used=set();branches=[]
    def walk(a,b):
        route=[a,b];used.add((min(a,b),max(a,b)))
        while degree[b]==2:
            row=sub.indices[sub.indptr[b]:sub.indptr[b+1]]
            c=int(row[0] if row[1]==a else row[1]);edge=(min(b,c),max(b,c))
            if edge in used:break
            used.add(edge);route.append(c);a,b=b,c
            if route[-1]==route[0]:break
        branches.append(ids[np.asarray(route)])
    for a in np.flatnonzero(degree!=2):
        for b in sub.indices[sub.indptr[a]:sub.indptr[a+1]]:
            if (min(a,int(b)),max(a,int(b))) not in used:walk(int(a),int(b))
    for a in range(len(ids)):
        for b in sub.indices[sub.indptr[a]:sub.indptr[a+1]]:
            if (min(a,int(b)),max(a,int(b))) not in used:walk(a,int(b))
    return branches


def prune_skeleton(points,graph,mask,minimum_length):
    active=mask.copy();history=[]
    for _ in range(len(mask)):
        branches=skeleton_branches(graph,active)
        degree=np.asarray((graph.astype(bool)@active.astype(np.int32))).ravel()
        remove=np.zeros(len(mask),dtype=bool)
        for b in branches:
            length=float(np.linalg.norm(np.diff(points[b],axis=0),axis=1).sum())
            if length>=minimum_length or b[0]==b[-1]:continue
            if degree[b[0]]==1 or degree[b[-1]]==1:
                remove[b[degree[b]<=2]]=True
        if not remove.any():return branches,history
        active[remove]=False;history.append(int(remove.sum()))
    raise RuntimeError('Skeleton pruning did not converge')


class PythonBackend:
    def __init__(self,directory,points,faces):
        self.directory=Path(directory);self.directory.mkdir(parents=True,exist_ok=True)
        centered=points-np.mean(points,axis=0)
        volume=np.einsum('ij,ij->',centered[faces[:,0]],
            np.cross(centered[faces[:,1]],centered[faces[:,2]]))/6
        if volume<=0:
            raise ValueError('Outward face winding is required for sulcal/gyral curvature; reverse inward faces before preprocessing')
        self.points=points;self.faces=faces;self.tools={};self.calls=[]
        self.graph,self.edges,self.lengths=weighted_graph(points,faces)
        self.indptr,self.links=vertex_links(faces,len(points))
        self.valid=np.ones(len(points),dtype=bool)

    def distance(self,mesh,mask,name):
        if not np.any(mask):return np.full(len(self.points),np.inf)
        return csgraph.dijkstra(self.graph,indices=np.flatnonzero(mask),min_only=True,directed=False)

    def skeleton(self,mask,method='distance',endpoint_percentile=5.):
        mask=mask&self.valid
        if not mask.any() or mask.all():
            return np.zeros(len(mask),dtype=bool),np.zeros(len(mask)),np.zeros(len(mask))
        if method=='hj-flux':
            depth,flux=hamilton_jacobi_flux(self.points,self.faces,self.graph,mask)
            threshold=float(np.percentile(flux[mask],endpoint_percentile))
            return thin_region_flux(mask,flux,threshold,self.indptr,self.links),depth,flux
        depth=self.distance(None,~mask,'boundary')
        return thin_region(mask,depth,self.indptr,self.links),depth,np.zeros(len(mask))


def segment(points,faces,atlas,feature,backend,config):
    from .preprocessing import adjacency,expand_ring,branch_map,write_branches,reconstruct_regions
    graph=adjacency(faces,len(points));ring2=expand_ring(graph)
    backend.valid=np.asarray(atlas)>0
    _,area,_=mesh_operators(points,faces)
    e=backend.edges
    if config.segmentation_energy=='robust-unit':
        costs,mixture=robust_unit_costs(feature,area,backend.valid)
        pairwise=np.full(len(e),config.graph_regularization)
    else:
        costs,mixture=gaussian_costs(feature,area,backend.valid)
        dual_width=(area[e[:,0]]+area[e[:,1]])/(3*backend.lengths)
        pairwise=config.graph_regularization*dual_width
        pairwise[~(backend.valid[e[:,0]]&backend.valid[e[:,1]])]=0
    initial,cut=potts_cut(costs,e,pairwise)
    labels=initial.copy();labels[~backend.valid]=0
    options=(config.skeleton_method,config.hj_endpoint_flux_percentile)
    cleaned=labels.copy()
    cleanup={'method':config.skeleton_relabeling,'flipped_gyral_vertices':0}
    if config.skeleton_relabeling=='pruned-branches':
        initial_gyral,_,_=backend.skeleton(labels==0,*options)
        retained_branches,initial_pruning=prune_skeleton(
            points,backend.graph,initial_gyral,config.gyral_pruning_mm)
        retained=np.zeros(len(points),dtype=bool)
        for branch in retained_branches:retained[branch]=True
        removed=initial_gyral&~retained
        if removed.any() and retained.any():
            removed_distance=backend.distance(None,removed,'removed_gyral')
            retained_distance=backend.distance(None,retained,'retained_gyral')
            flip=(labels==0)&backend.valid&(removed_distance<retained_distance)&(retained_distance>5)
            cleaned[flip]=1
            cleanup['flipped_gyral_vertices']=int(flip.sum())
        cleanup.update(initial_gyral_pruning=initial_pruning,
                       removed_gyral_skeleton_vertices=int(removed.sum()),
                       retained_gyral_skeleton_vertices=int(retained.sum()),
                       distance_rule='removed < retained and retained > 5 mm')
        cleaned_gyral,_,_=backend.skeleton(cleaned==0,*options)
        cleaned_sulcal,_,_=backend.skeleton(cleaned==1,*options)
        gb,_=prune_skeleton(points,backend.graph,cleaned_gyral,10.)
        sb,_=prune_skeleton(points,backend.graph,cleaned_sulcal,10.)
        labels=reconstruct_regions(cleaned,sb,gb,graph,valid=backend.valid,
            vertex_area=area,
            max_unsupported_area_mm2=config.reconstruction_max_unsupported_area_mm2)
        if not np.any(backend.valid&(labels==0)) or not np.any(backend.valid&(labels==1)):
            raise ValueError('Skeleton-guided region reconstruction removed a folding class')
        cleanup.update(reconstruction_gyral_branches=len(gb),reconstruction_sulcal_branches=len(sb),
                       reconstruction_max_unsupported_area_mm2=config.reconstruction_max_unsupported_area_mm2,
                       reconstruction_changed_vertices=int(np.sum(backend.valid&(labels!=cleaned))))
    sulcal,depth_s,flux=backend.skeleton(labels==1,*options)
    gyral,depth_g,_=backend.skeleton(labels==0,*options)
    branches,pruning=prune_skeleton(points,backend.graph,sulcal,config.sulcal_pruning_mm)
    gbranches,gpruning=prune_skeleton(points,backend.graph,gyral,config.gyral_pruning_mm)
    invalid=~backend.valid;excluded=invalid|np.asarray(ring2@invalid.astype(np.int32)).ravel().astype(bool)
    # Split at excluded vertices; never bridge across the medial-wall exclusion.
    final=[]
    for b in branches:
        kept=~excluded[b];runs=np.split(b,np.flatnonzero(kept[1:]!=kept[:-1])+1)
        final.extend(run for run in runs if len(run)>1 and not excluded[run[0]])
    branches=final
    history={'method':METHOD_VERSION,
        'feature':'minimum concave-positive principal curvature; negative of largest convex-positive principal curvature',
        'mixture':mixture,'graph_cut':cut,
        'segmentation_energy':config.segmentation_energy,
        'graph_regularization':float(config.graph_regularization),
        'skeleton_method':config.skeleton_method,
        'hj_endpoint_flux_percentile':config.hj_endpoint_flux_percentile if config.skeleton_method=='hj-flux' else None,
        'skeleton_relabeling':cleanup,
        'sulcal_deleted_vertices_per_pruning_pass':pruning,'gyral_deleted_vertices_per_pruning_pass':gpruning,
        'distance':'shortest paths along mesh edges in mm; not exact continuous geodesics',
        'skeleton':('one-ring average-flux priority, simple-point thinning and endpoint threshold'
                    if config.skeleton_method=='hj-flux' else 'distance-priority simple-point deletion on cyclic vertex links'),
        'legacy_flux_pruning_applied':False}
    (backend.directory/'python_methods.json').write_text(json.dumps(history,indent=2)+'\n')
    write_branches(backend.directory/'sulcal_branches.raw',branches)
    (backend.directory/'flux_pruning.json').write_text('[]\n')
    labels.astype('<i4').tofile(backend.directory/'segmentation.raw')
    depth=np.where(labels==1,depth_s,-depth_g);depth[~backend.valid]=0
    depth.astype('<f4').tofile(backend.directory/'signed_depth.raw')
    return dict(feature=feature,initial=initial,cleaned=cleaned,segmentation=labels,depth=depth,
        flux=np.asarray(flux,dtype=np.float32),sulcal_skeleton=branch_map(branches,len(points)),
        gyral_skeleton=branch_map(gbranches,len(points)),branches=branches,graph=graph,ring2=ring2,
        flux_history=[],method_details=history)
