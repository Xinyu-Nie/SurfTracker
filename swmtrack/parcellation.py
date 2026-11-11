"""Conservative endpoint-parcel reporting for a merged cortical target.

Tracking stops at a face boundary. A face gets a parcel only when all of its
gyral vertices agree; a face spanning parcels is reported as ambiguous (0).
This avoids moving endpoints or guessing from a nearby vertex.
"""
import csv
from pathlib import Path
import numpy as np


def cortical_labels(table):
    """Present positive annotation IDs, excluding named noncortical regions."""
    excluded={'unknown','???','corpuscallosum','medialwall','medial_wall'}
    labels=[int(item['id']) for item in table if item['vertices']>0 and item['id']>0
            and item['name'].lower().replace('-','').replace(' ','') not in excluded]
    if not labels:
        raise ValueError('Annotation contains no nonempty cortical parcels')
    if len(set(labels))!=len(labels):
        raise ValueError('Cortical annotation IDs must be distinct')
    return labels


def classify_endpoint_parcels(atlas,gyrus,faces,endpoint_faces):
    atlas=np.asarray(atlas,dtype=np.int64).ravel()
    gyrus=np.asarray(gyrus,dtype=bool).ravel()
    faces=np.asarray(faces,dtype=np.int64)
    ids=np.asarray(endpoint_faces,dtype=np.int64)
    if ids.ndim!=2 or ids.shape[1]!=2:
        raise ValueError('Endpoint faces must have shape (tracks, 2)')
    if len(atlas)!=len(gyrus) or np.any(ids<0) or np.any(ids>=len(faces)):
        raise ValueError('Endpoint faces or annotation dimensions are invalid')
    if not len(ids):return np.empty((0,2),np.int64),np.empty(0,np.uint8)
    triangles=faces[ids]
    values=np.where(gyrus[triangles],atlas[triangles],np.int64(0))
    upper=values.max(axis=-1)
    lower=np.where(values>0,values,np.iinfo(np.int64).max).min(axis=-1)
    if np.any(upper==0):
        raise ValueError('Endpoint face has no labeled gyral vertex')
    labels=np.where(lower==upper,upper,0).astype(np.int64)
    clear=(labels>0).all(axis=1)
    kind=np.full(len(ids),3,np.uint8)
    kind[clear&(labels[:,0]==labels[:,1])]=1
    kind[clear&(labels[:,0]!=labels[:,1])]=2
    return labels,kind


def write_endpoint_classification(directory,labels,kinds,lengths,trial_ids,table):
    """Write reports in tractogram order; lengths are integrated physical mm."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    pairs=np.sort(labels,axis=1)
    names={int(item['id']):item['name'] for item in table};names[0]='ambiguous'
    fields={'endpoint_a_id':labels[:,0],'endpoint_b_id':labels[:,1],
            'pair_low_id':pairs[:,0],'pair_high_id':pairs[:,1],
            'endpoint_class':kinds,'length_mm':np.asarray(lengths,dtype=np.float64),
            'accepted_trial_id':trial_ids}
    np.savez_compressed(directory/'endpoint_parcels.npz',**fields)
    with (directory/'endpoint_parcels.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['track_id','trial_id','length_mm','endpoint_a_id','endpoint_a_name',
            'endpoint_b_id','endpoint_b_name','pair_low_id','pair_high_id','endpoint_class'])
        for i,((a,b),(lo,hi),kind,length,trial) in enumerate(zip(labels,pairs,kinds,lengths,trial_ids)):
            w.writerow([i,int(trial),format(float(length),'.12g'),int(a),names[int(a)],int(b),names[int(b)],int(lo),int(hi),int(kind)])
    unique,counts=np.unique(pairs,axis=0,return_counts=True)
    with (directory/'parcel_pair_counts.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['parcel_low_id','parcel_low_name','parcel_high_id','parcel_high_name','tracks'])
        for (a,b),count in zip(unique,counts):w.writerow([int(a),names[int(a)],int(b),names[int(b)],int(count)])
    # A symmetric matrix, with each undirected connection counted once per
    # off-diagonal cell and same-parcel connections once on the diagonal.
    ids=[0]+sorted({int(item['id']) for item in table if item['id']>0 and item['vertices']>0})
    index={value:i for i,value in enumerate(ids)};matrix=np.zeros((len(ids),len(ids)),np.int64)
    for (a,b),count in zip(unique,counts):
        matrix[index[int(a)],index[int(b)]]=count
        if a!=b:matrix[index[int(b)],index[int(a)]]=count
    with (directory/'parcel_pair_matrix.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['id:name']+[f'{i}:{names[i]}' for i in ids])
        for i,row in zip(ids,matrix):w.writerow([f'{i}:{names[i]}',*map(int,row)])
    report={'tracks':len(labels),'same_parcel':int((kinds==1).sum()),
            'cross_parcel':int((kinds==2).sum()),'ambiguous_parcel':int((kinds==3).sum()),
            'endpoint_class_codes':{'1':'same_parcel','2':'cross_parcel','3':'ambiguous_parcel'},
            'assignment':'all gyral vertices of the terminal face must agree; otherwise parcel 0 (ambiguous)',
            'length_definition':'integrated trajectory length in physical mm; filtering performed during tracking',
            'exports':{name:str(directory/name) for name in ('endpoint_parcels.npz','endpoint_parcels.csv',
                'parcel_pair_counts.csv','parcel_pair_matrix.csv')}}
    return report,fields
