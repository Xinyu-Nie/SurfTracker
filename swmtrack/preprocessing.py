"""Sulcal/gyral preprocessing and the group's patch-selection workflow.

The default Python backend is an independent implementation of mesh methods.
The optional legacy-native backend requires separately supplied executables.
"""

from dataclasses import dataclass, asdict
from pathlib import Path
import json
import subprocess
import time
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from scipy.sparse.linalg import splu
from scipy.io import savemat
from numba import njit
from .io import read_surface, write_surface, sha256
from .geometry import vertex_normals
from .native import resolve_tools, check_native


@dataclass
class PreprocessConfig:
    """Explicit settings; no-config workflows use recommended().

    The bare constructor retains the historical cotan baseline for explicit
    comparisons. Use recommended() for the current workflow configuration.
    """
    backend: str = "python"
    curvature_method: str = "cotan"
    skeleton_method: str = "distance"
    skeleton_relabeling: str = "none"
    hj_endpoint_flux_percentile: float = 5.
    reconstruction_max_unsupported_area_mm2: float = 5.
    graph_regularization: float = .2
    branch_merge: str = "legacy"
    mesh_smoothing: float = .5
    feature_smoothing: float = .5
    smoothing_mass: str = "lumped"
    segmentation_energy: str = "area-weighted"
    quadric_radius_mm: float = 3.
    implicit_grid_spacing_mm: float = 1.
    implicit_derivative_sigma_mm: float = 1.
    implicit_gradient_fallback_threshold: float = .2
    implicit_band_mm: float = 4.
    implicit_candidate_faces: int = 64
    implicit_max_grid_voxels: int = 8_000_000
    deep_sulci_fraction: float = 0.
    gyral_pruning_mm: float = 25.
    sulcal_pruning_mm: float = 15.
    flux_count_threshold: float = -10.
    branch_include_fraction: float = .3
    patch_include_fraction: float = .3
    min_patch_branch_vertices: int = 15
    max_flux_pruning_iterations: int = 100
    between_gyri_seed_policy: str = "connected-intersection"
    within_gyrus_seed_policy: str = "intersection"

    @classmethod
    def recommended(cls, **overrides):
        """Validated workflow settings; explicit overrides take precedence."""
        return cls.agreement_candidate(**overrides)

    @classmethod
    def agreement_candidate(cls, **overrides):
        """Compatibility name for the settings evaluated against legacy labels.

        This retains existing Python HJ/cleanup and robust Gaussian fitting.
        Explicit keyword overrides take precedence over the candidate settings.
        """
        options = dict(curvature_method="implicit-grid", smoothing_mass="consistent",
                       implicit_grid_spacing_mm=1., implicit_derivative_sigma_mm=1.5,
                       implicit_band_mm=6., segmentation_energy="robust-unit",
                       graph_regularization=1., skeleton_method="hj-flux",
                       skeleton_relabeling="pruned-branches")
        options.update(overrides)
        return cls(**options)

    def validate(self):
        if self.backend not in ("python","legacy-native"):
            raise ValueError("backend must be python or legacy-native")
        if self.curvature_method not in ("cotan", "quadric", "implicit-grid"):
            raise ValueError("curvature_method must be cotan, quadric, or implicit-grid")
        if self.smoothing_mass not in ("lumped", "consistent"):
            raise ValueError("smoothing_mass must be lumped or consistent")
        if self.segmentation_energy not in ("area-weighted", "robust-unit"):
            raise ValueError("segmentation_energy must be area-weighted or robust-unit")
        if self.backend == "legacy-native" and (self.smoothing_mass != "lumped" or self.segmentation_energy != "area-weighted"):
            raise ValueError("smoothing_mass and segmentation_energy apply only to the Python backend")
        if self.skeleton_method not in ("distance", "hj-flux"):
            raise ValueError("skeleton_method must be distance or hj-flux")
        if self.skeleton_relabeling not in ("none", "pruned-branches"):
            raise ValueError("skeleton_relabeling must be none or pruned-branches")
        if not 0 <= self.hj_endpoint_flux_percentile <= 100:
            raise ValueError("hj_endpoint_flux_percentile must lie in [0, 100]")
        if self.reconstruction_max_unsupported_area_mm2 <= 0:
            raise ValueError("reconstruction_max_unsupported_area_mm2 must be positive")
        if self.backend == "legacy-native" and self.curvature_method != "cotan":
            raise ValueError("curvature_method applies only to the Python backend")
        if self.backend == "legacy-native" and (self.skeleton_method != "distance" or self.skeleton_relabeling != "none"):
            raise ValueError("skeleton_method and skeleton_relabeling apply only to the Python backend")
        if not np.isfinite([v for v in asdict(self).values() if isinstance(v,(float,int))]).all():
            raise ValueError("Preprocessing parameters must be finite")
        if self.branch_merge not in ("legacy","connected"):
            raise ValueError("branch_merge must be legacy or connected")
        if self.between_gyri_seed_policy not in ("connected-intersection", "intersection", "legacy-patches"):
            raise ValueError("between_gyri_seed_policy must be connected-intersection, intersection, or legacy-patches")
        if self.within_gyrus_seed_policy not in ("intersection", "legacy-patches"):
            raise ValueError("within_gyrus_seed_policy must be intersection or legacy-patches")
        if min(self.mesh_smoothing,self.feature_smoothing,self.graph_regularization) < 0:
            raise ValueError("Smoothing weights must be nonnegative")
        if min(self.quadric_radius_mm, self.implicit_grid_spacing_mm,
               self.implicit_derivative_sigma_mm, self.implicit_band_mm,
               self.implicit_gradient_fallback_threshold) <= 0:
            raise ValueError("Curvature radii, grid spacing, sigma, and band must be positive")
        if self.implicit_gradient_fallback_threshold > .5:
            raise ValueError("implicit_gradient_fallback_threshold must be at most 0.5")
        if self.implicit_band_mm < 3*self.implicit_derivative_sigma_mm:
            raise ValueError("Implicit band must cover at least three derivative sigmas")
        if self.implicit_candidate_faces < 1 or self.implicit_max_grid_voxels < 1:
            raise ValueError("Implicit candidate_faces and max_grid_voxels must be positive")
        if not 0 <= self.deep_sulci_fraction <= .3:
            raise ValueError("deep_sulci_fraction must lie in [0, .3]")
        if any(not 0 <= v <= 1 for v in (self.branch_include_fraction, self.patch_include_fraction)):
            raise ValueError("Inclusion fractions must lie in [0, 1]")
        if min(self.gyral_pruning_mm, self.sulcal_pruning_mm, self.min_patch_branch_vertices) < 0:
            raise ValueError("Pruning thresholds must be nonnegative")
        if self.max_flux_pruning_iterations < 0:
            raise ValueError("max_flux_pruning_iterations must be nonnegative")


def write_mni_obj(path, points, faces, *, legacy_rounding=False):
    """Internal native-tool interchange only; public geometry remains VTK/RAS mm."""
    points, faces = np.asarray(points), np.asarray(faces)
    tri = points[faces]
    normal = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
    normals = np.zeros_like(points, dtype=float)
    for corner in range(3):
        np.add.at(normals, faces[:,corner], normal)
    denominator = np.linalg.norm(normals, axis=1)+(1e-6 if legacy_rounding else 0)
    if np.any(denominator == 0):
        raise ValueError("Cannot write normals at isolated vertices")
    normals /= denominator[:,None]
    with Path(path).open("w") as stream:
        stream.write(f"P 0.3 0.7 0.5 100 1 {len(points)}\n")
        np.savetxt(stream, points, fmt="%.2f" if legacy_rounding else "%.17g")
        np.savetxt(stream, normals, fmt="%.2f" if legacy_rounding else "%.17g")
        stream.write(f"{len(faces)}\n0 1 1 1 1\n")
        np.savetxt(stream, (3*np.arange(1,len(faces)+1))[:,None], fmt="%d")
        np.savetxt(stream, faces, fmt="%d")


def read_raw(path, dtype, count=None):
    data = np.fromfile(path, dtype=dtype)
    if count is not None and len(data) != count:
        raise RuntimeError(f"{path}: expected {count} values, got {len(data)}")
    if not np.isfinite(data).all():
        raise RuntimeError(f"{path}: nonfinite values")
    return data


def read_branches(path, n_vertices):
    data = read_raw(path, "<i4")
    branches=[];offset=0
    while offset<len(data):
        count=int(data[offset]);offset+=1
        if count<0 or offset+count>len(data):
            raise ValueError("Malformed length-prefixed branch file")
        b=data[offset:offset+count].astype(np.int64);offset+=count
        if np.any(b<0) or np.any(b>=n_vertices):
            raise ValueError("Branch references a vertex outside the source mesh")
        if len(b):
            branches.append(b)
    return branches


def write_branches(path, branches):
    with Path(path).open("wb") as stream:
        for b in branches:
            np.r_[len(b),b].astype("<i4").tofile(stream)


def adjacency(faces, n):
    edges=np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    graph=sparse.coo_matrix((np.ones(2*len(edges),dtype=bool),
          (np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n)).tocsr()
    graph.sum_duplicates();graph.sort_indices()
    return graph


def expand_ring(graph):
    result=(graph@graph).astype(bool)
    result.setdiag(False);result.eliminate_zeros();result.sort_indices()
    return result


def components(graph, mask):
    ids=np.flatnonzero(mask);result=np.zeros(graph.shape[0],dtype=np.int32)
    if len(ids):
        _,labels=csgraph.connected_components(graph[ids][:,ids],directed=False)
        result[ids]=labels+1
    return result


def intersection_seed_domains(atlas, segmentation, graph, labels_a, labels_b):
    """Sulcal ROI intersection and components incident to the A/B interface.

    Connectivity uses only sulcal vertices within the two anatomical groups.
    A source is a sulcal vertex incident to an A--B mesh edge; the other end
    may be gyral. Every component touching a source is retained. Internal
    parcel edges within one group do not count as the A/B interface.
    """
    first=np.isin(np.asarray(atlas).ravel(),labels_a)
    second=np.isin(np.asarray(atlas).ravel(),labels_b)
    if np.any(first&second):
        raise ValueError("Between-gyrus label groups must be disjoint")
    raw=(np.asarray(segmentation).ravel()==1)&(first|second)
    component=components(graph,raw)
    row,col=graph.nonzero()
    crossing=(first[row]&second[col])|(second[row]&first[col])
    boundary=np.zeros(len(atlas),dtype=bool)
    boundary[row[crossing]]=True
    boundary&=raw
    retained=np.unique(component[boundary])
    connected=raw&np.isin(component,retained)
    sizes=np.bincount(component)
    details={"intersection_vertices":int(raw.sum()),
        "intersection_components":int(component.max(initial=0)),
        "atlas_interface_edges":int(crossing.sum()//2),
        "sulcal_boundary_vertices":int(boundary.sum()),
        "connected_vertices":int(connected.sum()),"retained_components":int(len(retained)),
        "components":[{"id":i,"vertices":int(sizes[i]),"retained":bool(i in retained)}
                      for i in range(1,len(sizes))]}
    return raw,connected,boundary,component,details


def label_gaussians(feature, fraction):
    """MATLAB quantile (Hazen), sample variance, and six-decimal text interface."""
    x=np.asarray(feature,dtype=np.float64)
    low,high=np.quantile(x,[fraction,(1+fraction)/2],method="hazen")
    groups=(x[(x>low)&(x<high)],x[x>high])
    parameters=[]
    for group in groups:
        if len(group)<2 or np.var(group,ddof=1)<=0:
            raise ValueError("Cannot estimate two nondegenerate feature classes")
        parameters.append([.5,float(group.mean()),float(group.var(ddof=1))])
    return np.array(parameters)


def branch_map(branches,n):
    result=np.zeros(n,dtype=np.int32)
    for i,b in enumerate(branches,1):
        result[b]=i
    return result


def ring_neighbors(graph,vertices):
    return graph[np.asarray(vertices)].indices


def connected_join(a,b,ring2):
    """Orient both branches so the join is a shared or two-ring neighbor pair."""
    for left,right in ((a,b),(a,b[::-1]),(a[::-1],b),(a[::-1],b[::-1])):
        if left[-1]==right[0] or right[0] in ring_neighbors(ring2,[left[-1]]):
            return np.r_[left,right]
    return None


def flux_prune(branches, flux, ring2, threshold=-10., max_iterations=100,merge_policy="legacy"):
    """Port the MATLAB post-pruning loop, including its removal-before-test order.

    Its branch score is mean(flux) * VERTEX COUNT, not physical branch length.
    Empty/no-leaf cases are guarded instead of indexing an empty MATLAB max.
    """
    branches=[b.copy() for b in branches];history=[]
    for iteration in range(max_iterations):
        mapping=branch_map(branches,len(flux));starts=[];ends=[];scores=[]
        for i,b in enumerate(branches,1):
            if len(b):
                starts.append(np.setdiff1d(np.unique(mapping[ring_neighbors(ring2,[b[0]])]),[0,i])-1)
                ends.append(np.setdiff1d(np.unique(mapping[ring_neighbors(ring2,[b[-1]])]),[0,i])-1)
                scores.append(float(np.mean(flux[b],dtype=np.float64)*len(b)))
            else:
                starts.append(np.array([],dtype=int));ends.append(np.array([],dtype=int));scores.append(0.)
        leaves=[i for i,b in enumerate(branches) if len(b) and bool(len(starts[i])) != bool(len(ends[i]))]
        if not leaves:
            history.append({"iteration":iteration,"stop":"no_leaf_guard"});break
        removed=leaves[int(np.argmax(np.asarray(scores)[leaves]))]
        history.append({"iteration":iteration,"removed_branch":removed,"score":scores[removed]})
        branches[removed]=np.array([],dtype=np.int64)
        for neighbors in (starts[removed],ends[removed]):
            if len(neighbors)==2:
                a,b=map(int,neighbors)
                if not len(branches[a]) or not len(branches[b]):
                    history.append({"iteration":iteration,"guard":"empty_merge_branch"});continue
                if merge_policy=="connected":
                    joined=connected_join(branches[a],branches[b],ring2)
                    if joined is None:
                        history.append({"iteration":iteration,"guard":"no_connected_join"});continue
                    branches[a]=joined
                elif merge_policy=="legacy":
                    if branches[a][-1] not in ring_neighbors(ring2,[branches[b][0]]):
                        branches[b]=branches[b][::-1]
                    branches[a]=np.r_[branches[a],branches[b]]
                else:raise ValueError("Unknown branch merge policy")
                branches[b]=np.array([],dtype=np.int64)
        if scores[removed]<threshold:
            break
    return [b for b in branches if len(b)],history


@njit(cache=True)
def _regular_subset(mask, faces, indptr, face_indices):
    """Legacy MeshSubset's ascending, sequential removal of pinched vertices."""
    result=mask.copy()
    changed=True
    while changed:
        changed=False
        for vertex in range(len(result)):
            if not result[vertex]:
                continue
            transitions=0
            for entry in range(indptr[vertex],indptr[vertex+1]):
                face=faces[face_indices[entry]]
                a,b=-1,-1
                for corner in range(3):
                    if face[corner]!=vertex:
                        if a<0: a=face[corner]
                        else: b=face[corner]
                if result[a]!=result[b]: transitions+=1
            if transitions>2:
                result[vertex]=False;changed=True
    used=np.zeros(len(mask),dtype=np.bool_)
    for face in faces:
        if result[face[0]] and result[face[1]] and result[face[2]]:
            used[face[0]]=True;used[face[1]]=True;used[face[2]]=True
    return used


def regular_subset(mask, faces):
    incidence=sparse.coo_matrix((np.ones(faces.size,dtype=bool),
              (faces.ravel(),np.repeat(np.arange(len(faces)),3))),shape=(len(mask),len(faces))).tocsr()
    return _regular_subset(np.asarray(mask,dtype=bool),np.asarray(faces,dtype=np.int64),incidence.indptr,incidence.indices)


def validate_full_surface(points,faces):
    """Native curvature/thinning kernels require a complete, closed hemisphere."""
    if points.ndim!=2 or points.shape[1]!=3 or not np.isfinite(points).all():
        raise ValueError("Preprocessing requires finite N x 3 world coordinates")
    if faces.ndim!=2 or faces.shape[1]!=3 or not len(faces) or faces.min()<0 or faces.max()>=len(points):
        raise ValueError("Invalid triangular surface")
    tri=points[faces]
    if np.any(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)<=1e-12):
        raise ValueError("Degenerate triangles must be repaired before preprocessing")
    if np.any(np.bincount(faces.ravel(),minlength=len(points))==0):
        raise ValueError("Isolated vertices must be removed before preprocessing")
    a=faces.ravel();b=np.roll(faces,-1,axis=1).ravel()
    low,high=np.minimum(a,b),np.maximum(a,b)
    order=np.lexsort((high,low))
    starts=np.flatnonzero(np.r_[True,(low[order[1:]]!=low[order[:-1]])|(high[order[1:]]!=high[order[:-1]]),True])
    if np.any(np.diff(starts)!=2):
        raise ValueError("Preprocessing requires a closed manifold hemisphere; apply tracking ROIs afterwards")
    left,right=order[starts[:-1]],order[starts[:-1]+1]
    if np.any(a[left]!=b[right]) or np.any(b[left]!=a[right]):
        raise ValueError("Surface face winding is inconsistent")


class LegacyNative:
    def __init__(self, directory, overrides=None):
        self.directory=Path(directory).resolve();self.directory.mkdir(parents=True,exist_ok=True)
        self.tools=resolve_tools(overrides)
        self.overrides=overrides
        self.checked=False
        history=self.directory/"native_calls.json"
        self.calls=json.loads(history.read_text()) if history.exists() else []

    def run(self, tool, *arguments):
        if not self.checked:
            check_native(self.overrides)
            self.checked=True
        executable=self.tools[tool]
        if not executable.is_file():
            raise FileNotFoundError(f"Missing legacy-native kernel {tool}: {executable}; supply tool overrides")
        log=self.directory/f"{len(self.calls):02d}_{tool}.log"
        command=[str(executable),*[str(v) for v in arguments]]
        started=time.perf_counter()
        with log.open("w") as stream:
            subprocess.run(command,cwd=self.directory,stdout=stream,stderr=subprocess.STDOUT,check=True)
        self.calls.append({"tool":tool,"command":command,"seconds":time.perf_counter()-started,"log":str(log)})
        (self.directory/"native_calls.json").write_text(json.dumps(self.calls,indent=2)+"\n")

    def distance(self, mesh, mask, name):
        mask=np.asarray(mask,dtype=bool)
        if not mask.any():
            return np.full(len(mask),np.inf)
        source=self.directory/f"{name}.raw";destination=self.directory/f"{name}_distance.raw"
        mask.astype("<i4").tofile(source)
        self.run("GeoDist",mesh,source,destination)
        return read_raw(destination,"<f4",len(mask)).astype(float)

    def skeleton(self, mesh, segmentation, region, prefix):
        source=self.directory/f"{prefix}_seg.raw";destination=self.directory/f"{prefix}.raw"
        np.asarray(segmentation,dtype="<i4").tofile(source)
        selected=(np.asarray(segmentation)==(0 if region==1 else 1))
        if not selected.any() or selected.all():
            return np.zeros(len(selected),dtype=np.int32),np.zeros(len(selected),dtype=np.float32)
        self.run("CortexSkeleton",mesh,source,region,.1,destination)
        values=read_raw(destination,"<i4",len(selected))
        flux=read_raw(self.directory/"flux.raw","<f4",len(selected))
        flux.tofile(self.directory/f"{prefix}_flux.raw")
        return values,flux

    def prune(self, mesh, skeleton, threshold, prefix):
        if not np.any(skeleton):
            return []
        source=self.directory/f"{prefix}_skeleton.raw";graph=self.directory/f"{prefix}_graph.raw";pruned=self.directory/f"{prefix}_pruned.raw"
        np.asarray(skeleton,dtype="<i4").tofile(source)
        self.run("SkeletonBranchExtraction",mesh,source,graph)
        self.run("SkeletonPruning",mesh,graph,threshold,0,pruned)
        return read_branches(pruned,len(skeleton))


def feature_k1(points,faces,backend,config):
    """Original implicit curvature kernel plus Python Tikhonov sparse solves."""
    mesh=backend.directory/"input.obj";write_mni_obj(mesh,points,faces)
    qfile,kfile=backend.directory/"Q.raw",backend.directory/"K.raw"
    backend.run("MeshLaplaceBeltramiEigen",mesh,qfile,kfile)
    matrices=[]
    for path in (qfile,kfile):
        data=read_raw(path,"<f8").reshape(-1,3)
        matrices.append(sparse.coo_matrix((data[:,2],(data[:,0].astype(int),data[:,1].astype(int))),shape=(len(points),len(points))).tocsc())
    q,k=matrices
    smoothed=splu(k+config.mesh_smoothing*q).solve(k@points)
    smooth_obj=backend.directory/"smoothed_legacy_rounded.obj"
    write_mni_obj(smooth_obj,smoothed,faces,legacy_rounding=True)
    grid=backend.directory/"smoothed.dtg"
    backend.run("Mesh2DTGrid",smooth_obj,1,4,grid)
    curvature,mc,gc=[backend.directory/name for name in ("curvature.raw","mean_curvature.raw","gaussian_curvature.raw")]
    backend.run("ImplicitCurvatureEstimation",grid,smooth_obj,curvature,mc,gc,0)
    mean=read_raw(mc,"<f4",len(points)).astype(float);gauss=read_raw(gc,"<f4",len(points)).astype(float)
    raw=mean+np.sqrt(np.maximum(mean*mean-gauss,0))
    result=splu(k+config.feature_smoothing*q).solve(k@raw).astype(np.float32)
    result.tofile(backend.directory/"k1.raw")
    np.savez_compressed(backend.directory/"feature_intermediates.npz",smoothed=smoothed,mean_curvature=mean,gaussian_curvature=gauss,k1_unfiltered=raw,k1=result)
    return result


def reconstruct_regions(segmentation,sulcal_branches,gyral_branches,graph,valid=None,
                        vertex_area=None,max_unsupported_area_mm2=None):
    """Port ReconSulcalGyralRegionFromPruning's component-based region growing."""
    original=np.asarray(segmentation,dtype=np.int32)
    labels=original.copy()
    sulcal=np.concatenate(sulcal_branches) if sulcal_branches else np.array([],dtype=int)
    gyral=np.concatenate(gyral_branches) if gyral_branches else np.array([],dtype=int)
    if valid is not None:
        valid=np.asarray(valid,dtype=bool)
        sulcal=sulcal[valid[sulcal]]
        gyral=gyral[valid[gyral]]
    labels[sulcal]=3;labels[gyral]=2
    # Seeds themselves have labels 3/2; grow only into original class 1/0.
    for seed,class_value,grown_value in ((sulcal,1,3),(gyral,0,2)):
        mask=(labels==class_value)|(labels==grown_value)
        if valid is not None:mask &= valid
        cc=components(graph,mask);touched=np.unique(cc[seed]);touched=touched[touched>0]
        labels[np.isin(cc,touched)&(labels==class_value)]=grown_value
    labels[labels==0]=3;labels[labels==1]=2
    labels-=2
    if valid is not None:labels[~valid]=0
    if max_unsupported_area_mm2 is not None:
        if valid is None or vertex_area is None:
            raise ValueError('Area-limited reconstruction needs valid mask and vertex areas')
        for value in (0,1):
            initial_mask=valid&(original==value)
            component=components(graph,initial_mask)
            size=np.bincount(component,weights=np.where(initial_mask,vertex_area,0),
                             minlength=int(component.max())+1)
            preserve=(component>0)&(size[component]>max_unsupported_area_mm2)&initial_mask
            labels[preserve]=value
    return labels


def segment_sulci(points,faces,atlas,feature,backend,config):
    n=len(points);graph=adjacency(faces,n);ring2=expand_ring(graph)
    mesh=backend.directory/"input.obj"
    if not mesh.exists(): write_mni_obj(mesh,points,faces)
    feature_path=backend.directory/"k1.raw";np.asarray(feature,dtype="<f4").tofile(feature_path)
    parameters=label_gaussians(feature,config.deep_sulci_fraction)
    parameter_path=backend.directory/"k1_gaussians.txt"
    parameter_path.write_text("".join("1\n"+"\n".join(f"{v:.6f}" for v in row)+"\n" for row in parameters))
    initial_path=backend.directory/"initial_segmentation.raw"
    backend.run("MeshSegmentation",mesh,feature_path,parameter_path,1,1,initial_path)
    initial=read_raw(initial_path,"<i4",n)
    gyral,_=backend.skeleton(mesh,initial,1,"initial_gyral")
    gyral_branches=backend.prune(mesh,gyral,config.gyral_pruning_mm,"initial_gyral")
    retained=branch_map(gyral_branches,n)>0
    removed=(~retained)&(gyral>0)
    da=backend.distance(mesh,removed,"removed_gyral")
    db=backend.distance(mesh,retained,"retained_gyral")
    cleaned=initial.copy();cleaned[(initial==0)&(da<db)&(db>5)]=1
    gyral,_=backend.skeleton(mesh,cleaned,1,"cleaned_gyral")
    sulcal,_=backend.skeleton(mesh,cleaned,0,"cleaned_sulcal")
    gb=backend.prune(mesh,gyral,10,"cleaned_gyral")
    sb=backend.prune(mesh,sulcal,10,"cleaned_sulcal")
    segmentation=reconstruct_regions(cleaned,sb,gb,graph)
    if np.all(segmentation==0) or np.all(segmentation==1):
        raise ValueError("Reconstructed sulcal/gyral segmentation has only one class; review features and pruning thresholds")
    d_sulcal=backend.distance(mesh,segmentation==1,"sulcal_region")
    d_gyral=backend.distance(mesh,segmentation==0,"gyral_region")
    depth=np.where(segmentation==0,-d_sulcal,d_gyral)
    sulcal,flux=backend.skeleton(mesh,segmentation,0,"final_sulcal")
    invalid=np.asarray(atlas)==0
    excluded=invalid|np.asarray(ring2@invalid.astype(np.int32)).ravel().astype(bool)
    sulcal[excluded]=0
    branches=backend.prune(mesh,sulcal,config.sulcal_pruning_mm,"final_sulcal")
    before=[b.copy() for b in branches]
    branches,history=flux_prune(branches,flux,ring2,config.flux_count_threshold,config.max_flux_pruning_iterations,config.branch_merge)
    segmentation.astype("<i4").tofile(backend.directory/"segmentation.raw")
    depth.astype("<f4").tofile(backend.directory/"signed_depth.raw")
    write_branches(backend.directory/"sulcal_branches.raw",branches)
    cells=np.empty((1,len(branches)),dtype=object)
    for i,b in enumerate(branches):cells[0,i]=b[:,None]
    savemat(backend.directory/"sulcal_branches.mat",{"BranchSet":cells})
    (backend.directory/"flux_pruning.json").write_text(json.dumps(history,indent=2)+"\n")
    return {"feature":feature,"initial":initial,"cleaned":cleaned,"segmentation":segmentation,"depth":depth,
            "flux":flux,"sulcal_skeleton":sulcal,"gyral_skeleton":gyral,"branches":branches,
            "pre_flux_branches":before,"graph":graph,"ring2":ring2,"gaussians":parameters,"flux_history":history}


def make_patches(points,faces,atlas,segmentation,branches,graph,ring2,backend,config,labels_a,labels_b=None):
    """Port one-/two-gyrus branch grouping, geodesic competition and patch tests."""
    first=np.isin(atlas,labels_a);second=np.isin(atlas,labels_b) if labels_b is not None else None
    mapping=branch_map(branches,len(points));ring4=expand_ring(ring2)
    selected=[]
    for i,b in enumerate(branches):
        count=int(first[b].sum())+(int(second[b].sum()) if second is not None else 0)
        if count>len(b)*config.branch_include_fraction:selected.append(i)
    # Keep only selected branch-to-branch edges; match the four-ring convention.
    row=[];column=[]
    for i in selected:
        neighbors=np.setdiff1d(np.unique(mapping[ring_neighbors(ring4,branches[i])]),[0,i+1])-1
        for j in np.intersect1d(neighbors,selected):row.append(i);column.append(int(j))
    branch_graph=sparse.coo_matrix((np.ones(len(row),dtype=bool),(row,column)),shape=(len(branches),len(branches))).tocsr()
    picked=np.zeros(len(branches),dtype=bool);picked[selected]=True
    branch_cc=components(branch_graph,picked)
    sulcal_cc=components(graph,segmentation==1)
    patches=[];decisions=[]
    for component in range(1,int(branch_cc.max(initial=0))+1):
        ids=np.flatnonzero(branch_cc==component)
        vertices=np.concatenate([branches[i] for i in ids])
        if len(vertices)<config.min_patch_branch_vertices:continue
        labels=np.unique(sulcal_cc[vertices])
        if len(labels)!=1 or labels[0]==0:
            decisions.append({"component":component,"skip":"skeleton spans multiple or non-sulcal components","labels":labels.tolist()});continue
        sulcus=sulcal_cc==labels[0]
        a=np.zeros(len(points),dtype=bool);a[vertices]=True
        b=(mapping>0)&sulcus;b[vertices]=False
        if b.any():
            da=backend.distance(backend.directory/"input.obj",a,f"patch_{component}_selected")
            db=backend.distance(backend.directory/"input.obj",b,f"patch_{component}_other")
            mask=sulcus&(da<db)
        else:mask=sulcus
        mask=mask|np.asarray(ring2@mask.astype(np.int32)).ravel().astype(bool)
        mask=regular_subset(mask,faces)
        cc=components(graph,mask)
        sizes=np.bincount(cc);sizes[0]=0
        if not sizes.any():continue
        mask=cc==int(np.argmax(sizes))
        fraction_a=float(first[mask].mean());fraction_b=float(second[mask].mean()) if second is not None else None
        accepted=fraction_a>config.patch_include_fraction and (second is None or fraction_b>config.patch_include_fraction)
        decisions.append({"component":component,"branch_ids":ids.tolist(),"vertices":int(mask.sum()),
                          "fraction_a":fraction_a,"fraction_b":fraction_b,"accepted":bool(accepted)})
        if accepted:patches.append(mask)
    return patches,decisions


def preprocess(surface,output_dir,*,config=None,labels_a=(14423100,),labels_b=(1316060,),
               atlas_array="aparc",tool_overrides=None,input_feature=None,stop_after="patches",reuse_segmentation=None,
               whole_hemisphere=False):
    config=config or PreprocessConfig.recommended();config.validate()
    if whole_hemisphere and (labels_b is not None or config.within_gyrus_seed_policy!="intersection"):
        raise ValueError("Whole-hemisphere preprocessing requires one merged target and intersection seeds")
    started=time.perf_counter();output=Path(output_dir).resolve();output.mkdir(parents=True,exist_ok=True)
    points,faces,data=read_surface(surface)
    validate_full_surface(points,faces)
    if atlas_array not in data:raise ValueError(f"Missing atlas point array {atlas_array}")
    atlas=np.asarray(data[atlas_array]).ravel()
    if len(atlas)!=len(points):raise ValueError("Atlas length differs from surface")
    if not np.isfinite(atlas).all() or np.any(atlas<0) or np.any(atlas!=np.floor(atlas)):
        raise ValueError("Atlas must contain nonnegative integer IDs, with 0 for unknown vertices")
    if stop_after=="patches":
        if not len(labels_a):
            raise ValueError("At least one anatomical label is required")
        requested=list(labels_a)+(list(labels_b) if labels_b is not None else [])
        if any(label<=0 or not np.any(atlas==label) for label in requested):
            raise ValueError("Each requested anatomical label must be positive and present in the atlas")
        if len(set(labels_a))!=len(labels_a) or (labels_b is not None and len(set(labels_b))!=len(labels_b)):
            raise ValueError("Anatomical labels within each group must be distinct")
        if labels_b is not None and (not len(labels_b) or np.intersect1d(labels_a,labels_b).size):
            raise ValueError("Between-gyrus label groups must be nonempty and disjoint")
    if config.backend=="python":
        if tool_overrides:raise ValueError("Tool overrides apply only to legacy-native preprocessing")
        from .open_preprocessing import PythonBackend,curvature_feature,curvature_feature_options,segment as python_segment,METHOD_VERSION
        backend=PythonBackend(output/"intermediates",points,faces)
    else:
        backend=LegacyNative(output/"intermediates",tool_overrides)
        mesh=backend.directory/"input.obj";write_mni_obj(mesh,points,faces)
    cached=None
    if reuse_segmentation is not None:
        if input_feature is not None:raise ValueError("Choose an input feature or a segmentation cache, not both")
        cache=Path(reuse_segmentation).resolve()
        if cache==output:raise ValueError("Write new patch selections to a separate output directory")
        parent=json.loads((cache/"preprocessing.json").read_text())
        if config.backend=="python" and parent.get("method_version")!=METHOD_VERSION:
            raise ValueError("Segmentation cache method version differs; recompute preprocessing with the current Python backend")
        if parent["surface_sha256"]!=sha256(surface) or parent["atlas_array"]!=atlas_array:
            raise ValueError("Segmentation cache belongs to a different surface or atlas")
        patch_fields={"branch_include_fraction","patch_include_fraction","min_patch_branch_vertices","branch_merge",
                      "between_gyri_seed_policy","within_gyrus_seed_policy"}
        previous_config=asdict(PreprocessConfig(backend="legacy-native"));previous_config.update(parent["config"])
        if any(previous_config[key]!=value for key,value in asdict(config).items() if key not in patch_fields):
            raise ValueError("Segmentation parameters differ from the supplied cache")
        with np.load(cache/"segmentation.npz",allow_pickle=False) as stored:
            cached={name:stored[name].copy() for name in stored.files}
        if any(len(v)!=len(points) or not np.isfinite(v).all() for v in cached.values()):
            raise ValueError("Invalid per-vertex segmentation cache arrays")
        cached["branches"]=read_branches(cache/"intermediates/sulcal_branches.raw",len(points))
        cached["graph"]=adjacency(faces,len(points));cached["ring2"]=expand_ring(cached["graph"])
        cached["flux_history"]=json.loads((cache/"intermediates/flux_pruning.json").read_text())
        before=cache/"pre_flux_branches.raw"
        if not before.exists():before=cache/"intermediates/final_sulcal_pruned.raw"
        if before.exists() and config.backend=="legacy-native":
            cached["pre_flux_branches"]=read_branches(before,len(points))
            cached["branches"],cached["flux_history"]=flux_prune(cached["pre_flux_branches"],cached["flux"],cached["ring2"],
                config.flux_count_threshold,config.max_flux_pruning_iterations,config.branch_merge)
        elif parent["config"].get("branch_merge","legacy")!=config.branch_merge:
            raise ValueError("This cache lacks pre-flux branch data; reuse the original segmentation cache")
        feature=cached["feature"]
        write_branches(backend.directory/"sulcal_branches.raw",cached["branches"])
        (backend.directory/"flux_pruning.json").write_text(json.dumps(cached["flux_history"],indent=2)+"\n")
        np.asarray(feature,dtype="<f4").tofile(backend.directory/"k1.raw")
    else:
        if input_feature is not None:feature=read_raw(input_feature,"<f4",len(points))
        elif config.backend=="python":
            if config.curvature_method == "cotan":
                feature,intermediates=curvature_feature(points,faces,config.mesh_smoothing,config.feature_smoothing,
                                                      smoothing_mass=config.smoothing_mass)
            else:
                feature,intermediates,estimator_details=curvature_feature_options(points,faces,config)
                (backend.directory/"curvature_estimator.json").write_text(json.dumps(estimator_details,indent=2)+"\n")
            np.savez_compressed(backend.directory/"feature_intermediates.npz",**intermediates)
        else:feature=feature_k1(points,faces,backend,config)
        np.asarray(feature,dtype="<f4").tofile(backend.directory/"k1.raw")
    arrays={**data,"k1":feature};result={"feature":feature}
    metadata={"surface":str(Path(surface).resolve()),"surface_sha256":sha256(surface),"backend":config.backend,
              "coordinate_system":"RAS+ mm","config":asdict(config),"atlas_array":atlas_array,
              "labels_a":list(labels_a),"labels_b":list(labels_b) if labels_b is not None else None,
              "scope":"MATLAB-free orchestration with packaged original native kernels",
              "curvature_geometry":"Tikhonov-smoothed and rounded to 0.01 mm only for legacy implicit curvature; tracking surface remains unchanged",
              "endpoint_definition":"annotation region intersected with reconstructed gyral class, excluding accepted seed patches",
              "legacy_flux_score":"mean flux times branch vertex count; remove-before-threshold-test order retained"}
    if input_feature is not None:
        metadata["input_feature"]={"path":str(Path(input_feature).resolve()),"sha256":sha256(input_feature)}
    if reuse_segmentation is not None:
        metadata["reused_segmentation"]={"directory":str(cache),"parent_metadata_sha256":sha256(cache/"preprocessing.json")}
    if stop_after!="features":
        result=cached if cached is not None else (python_segment(points,faces,atlas,feature,backend,config)
            if config.backend=="python" else segment_sulci(points,faces,atlas,feature,backend,config))
        if "pre_flux_branches" in result:
            write_branches(output/"pre_flux_branches.raw",result["pre_flux_branches"])
        arrays.update({"sulcal_segmentation":result["segmentation"],"signed_sulcal_depth_mm":result["depth"],
                       "sulcal_skeleton":result["sulcal_skeleton"],"gyral_skeleton":result["gyral_skeleton"],
                       "sulcal_flux":result["flux"],"sulcal_branch_id":branch_map(result["branches"],len(points))})
        if config.backend=="python" and config.skeleton_method!="hj-flux":arrays.pop("sulcal_flux",None)
        np.savez_compressed(output/"segmentation.npz",**{k:result[k] for k in ("feature","initial","cleaned","segmentation","depth","flux","sulcal_skeleton","gyral_skeleton")})
        if stop_after=="patches":
            # These arrays belong to this selection; do not carry an earlier
            # within/between selection's masks into a newly selected region.
            for name in ("seed","gyrus","gyrus_a","gyrus_b","roi","seed_intersection",
                         "seed_connected","seed_component_id","interface_sulcal_vertices"):
                arrays.pop(name,None);arrays.pop(name+"_RGB",None)
            policy=config.between_gyri_seed_policy if labels_b is not None else config.within_gyrus_seed_policy
            if policy=="legacy-patches":
                patches,decisions=make_patches(points,faces,atlas,result["segmentation"],result["branches"],
                    result["graph"],result["ring2"],backend,config,labels_a,labels_b)
                seed=np.any(patches,axis=0) if patches else np.zeros(len(points),dtype=bool)
                raw=seed
            elif labels_b is None:
                region=np.isin(atlas,labels_a)
                if whole_hemisphere and "cortex" in data:
                    region &= np.asarray(data["cortex"]).ravel()>0
                if not np.isin(result["segmentation"][region],[0,1]).all():
                    raise ValueError("Selected ROI contains an invalid sulcal/gyral segmentation class")
                raw=(result["segmentation"]==1)&region
                seed=raw
                component=components(result["graph"],raw)
                sizes=np.bincount(component)
                decisions=[{"id":i,"vertices":int(sizes[i]),"accepted":True}
                           for i in range(1,len(sizes))]
                patches=[component==i for i in range(1,len(sizes))]
                arrays.update(seed_intersection=raw,seed_component_id=component)
                metadata.update(seed_selection={"intersection_vertices":int(raw.sum()),
                    "intersection_components":len(decisions),"retained_components":len(decisions),
                    "components":decisions},seed_selection_version="sulcal-roi-intersection-v1")
            else:
                raw,connected,boundary,component,details=intersection_seed_domains(
                    atlas,result["segmentation"],result["graph"],labels_a,labels_b)
                seed=connected if policy=="connected-intersection" else raw
                patches=[(component==i)&seed for i in np.unique(component[seed])]
                decisions=details["components"]
                arrays.update(seed_intersection=raw,seed_connected=connected,
                    seed_component_id=component,interface_sulcal_vertices=boundary)
                metadata.update(seed_selection=details,seed_selection_version="sulcal-roi-interface-v1")
            gyral=(result["segmentation"]==0)&~seed
            a=gyral&np.isin(atlas,labels_a)
            if whole_hemisphere:
                a &= region
            b=gyral&np.isin(atlas,labels_b) if labels_b is not None else np.zeros(len(points),dtype=bool)
            # Keep the full sulcal intersection traversable in BOTH modes.
            # Changing seed eligibility must not also change geometry or targets.
            roi=raw|a|b
            arrays.update({"seed":seed,"gyrus_a":a,"gyrus_b":b,"roi":roi})
            if labels_b is None:arrays["gyrus"]=a
            for i,mask in enumerate(patches,1):
                kept=np.all(mask[faces],axis=1);ids=np.unique(faces[kept]);mapping=np.full(len(points),-1,dtype=int);mapping[ids]=np.arange(len(ids))
                if not len(ids):continue
                write_surface(output/f"patch_{i:03d}.vtk",points[ids],mapping[faces[kept]],{k:v[ids] for k,v in arrays.items()})
                np.save(output/f"patch_{i:03d}_source_vertex_ids.npy",np.asarray(data.get("source_vertex_id",np.arange(len(points))))[ids])
            metadata.update({"patch_count":len(patches),"patch_decisions":decisions,"seed_vertices":int(seed.sum()),
                             "seed_policy":policy,"gyrus_a_vertices":int(a.sum()),"gyrus_b_vertices":int(b.sum()),"roi_vertices":int(roi.sum()),
                             "endpoint_policy":"within_gyrus" if labels_b is None else "between_gyri",
                             "roi_definition":"full sulcal intersection and gyral endpoints" if policy!="legacy-patches" else "accepted patches and gyral endpoints",
                             "tracking_ready":bool(seed.any() and a.any() and (labels_b is None or b.any()))})
            if labels_b is None:
                metadata["endpoint_definition"]="one binary union of the selected annotation regions intersected with gyral class 0, excluding seed vertices"
            if whole_hemisphere:
                metadata.update(tracking_scope="whole-hemisphere",cortex_mask_used="cortex" in data,
                    roi_definition="all selected cortical annotation vertices; medial wall excluded")
    write_surface(output/"preprocessed.world.vtk",points,faces,arrays)
    metadata["native_tools"]={name:{"path":str(path),"sha256":sha256(path)} for name,path in backend.tools.items() if any(c["tool"]==name for c in backend.calls)}
    metadata["native_calls"]=backend.calls;metadata["elapsed_seconds"]=time.perf_counter()-started
    if config.backend=="python":
        metadata.update(scope="Independent mesh preprocessing using NumPy/SciPy/Numba; no legacy executable calls",
            method_version=METHOD_VERSION,curvature_geometry={
                "cotan":"Cotangent mean curvature and angle-defect Gaussian curvature on smoothed mesh",
                "quadric":"Geodesic-neighborhood quadric fit on smoothed mesh",
                "implicit-grid":"Signed-distance grid and tangent-plane shape operator from smoothed mesh"}[
                    config.curvature_method]+"; tracking geometry unchanged",
            feature_definition="smoothed minimum concave-positive principal curvature H-sqrt(max(H^2-K,0)); negative of largest convex-positive principal curvature",
            legacy_flux_score=None,legacy_only_parameters=["branch_merge","deep_sulci_fraction","flux_count_threshold","max_flux_pruning_iterations"],
            distances="Mesh-edge shortest paths in mm",skeleton=(
                "One-ring average-flux priority and simple-point thinning with physical leaf pruning"
                if config.skeleton_method=="hj-flux" else "Distance-priority vertex-link thinning with physical leaf pruning"),
            source_sha256_python_backend=sha256(Path(__file__).with_name("open_preprocessing.py")))
        if config.curvature_method != "cotan":
            metadata["source_sha256_curvature_estimators"]=sha256(Path(__file__).with_name("curvature_estimators.py"))
        if "method_details" in result:metadata["method_details"]=result["method_details"]
    metadata["source_sha256"]=sha256(Path(__file__))
    (output/"preprocessing.json").write_text(json.dumps(metadata,indent=2)+"\n")
    return metadata
