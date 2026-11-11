"""Optional ANTsPy subprocess worker; imports no core swmtrack dependencies."""
import argparse
import json
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--probe',action='store_true');parser.add_argument('--request')
    args=parser.parse_args()
    import ants
    if not hasattr(ants,'config') or not hasattr(ants.config,'set_ants_deterministic'):
        raise RuntimeError('ANTsPy >= 0.6.3 is required for explicit registration seed control')
    if args.probe:
        print(json.dumps({'python':str(Path(sys.executable).absolute()),'ants_version':ants.__version__}));return
    import numpy as np
    import pandas as pd
    request=json.loads(Path(args.request).read_text());out=Path(request['output'])
    ants.config.set_ants_deterministic(on=False,seed_value=request['random_seed'])
    fixed=ants.image_read(str(out/'mean_b0.nii.gz'));moving=ants.image_read(str(out/'t1_brain.nii.gz'))
    fixed_mask=ants.image_read(str(out/'fixed_mask.nii.gz')) if request['provided_mask'] else ants.get_mask(fixed)
    if not np.any(fixed_mask.numpy()>0):raise ValueError('Automatic fixed-image mask is empty; supply --dwi-mask')
    ants.image_write(fixed_mask,str(out/'fixed_mask.nii.gz'))
    moving_mask=moving>0;ants.image_write(moving_mask,str(out/'moving_mask.nii.gz'))
    result=ants.registration(fixed=fixed,moving=moving,mask=fixed_mask,moving_mask=moving_mask,
        type_of_transform='Rigid',initial_transform='Identity' if request['initialization']=='header' else None,
        outprefix=str(out/'t1_to_b0_'),aff_metric='mattes',aff_sampling=32,aff_random_sampling_rate=.25,
        aff_iterations=(1000,500,250,100),aff_shrink_factors=(8,4,2,1),aff_smoothing_sigmas=(3,2,1,0),
        smoothing_in_mm=True,mask_all_stages=True,singleprecision=False,verbose=False)
    transforms=result['fwdtransforms']
    if len(transforms)!=1 or not transforms[0].endswith('.mat'):raise ValueError('Expected one rigid ANTs affine')
    native=ants.read_transform(transforms[0],precision='double')
    a=native.parameters[:9].reshape(3,3);t=native.parameters[9:];center=native.fixed_parameters
    pull=np.eye(4);pull[:3,:3]=a;pull[:3,3]=t+center-a@center
    flip=np.diag([-1.,-1.,1.,1.]);matrix=flip@np.linalg.inv(pull)@flip
    points=np.array([[-40.,-20.,30.],[-10.,20.,50.],[30.,-60.,10.],[0.,0.,0.]])
    actual=ants.apply_transforms_to_points(3,pd.DataFrame(points*[-1.,-1.,1.],columns=['x','y','z']),
        transforms,whichtoinvert=[True])[['x','y','z']].to_numpy()*[-1.,-1.,1.]
    expected=points@matrix[:3,:3].T+matrix[:3,3]
    error=float(np.max(np.abs(actual-expected)))
    if error>3e-5:raise ValueError('ANTs LPS/RAS point-transform convention check failed')
    np.savetxt(out/'scanner_ras_to_fod_world.txt',matrix,fmt='%.17g')
    ants.image_write(result['warpedmovout'],str(out/'t1_in_b0.nii.gz'))
    report={**request,'ants_version':ants.__version__,'python':str(Path(sys.executable).absolute()),
        'fixed_mask_source':'provided' if request['provided_mask'] else 'ANTs get_mask intensity mask',
        'scanner_ras_to_fod_world':matrix.tolist(),'point_application_max_error_mm':error,
        'point_mapping':'inverse of fixed-to-moving image pull map, LPS conjugated to RAS',
        'fod_resampled':False,'anatomical_alignment_verified':False}
    (out/'registration.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
