function WriteTrackVisFile(LineInd,PtVec,dim, voxel_size,fname)
%WriteTrackVisFile(LineInd,PtVec,dim, voxel_size,fname)
%
%PtVec: points, assume in voxel xyz coordinate
% dim: 1x3 array, size of the image array
% voxel_size: 1x3 array, voxel scale
%
%Yonggang Shi
% USC Stevens Neuroimaging and Informatics Institute
% Laboratory of Neuro Imaging (LONI)
% Keck School of Medicine of USC
% University of Southern California


fid = fopen(fname,'w');
id_string = [84  82   65    67    75      0 ];
fwrite(fid,id_string,'char');

fwrite(fid,dim,'int16');
fwrite(fid,voxel_size,'float32');

origin = [0 0 0];
fwrite(fid,origin,'float32');
n_s = 0;
fwrite(fid,n_s,'int16');
s_name = zeros(200,1);
fwrite(fid,s_name,'char');
n_p = 0;
fwrite(fid,n_p,'int16');
p_name = zeros(200,1);
fwrite(fid,p_name,'char');
vox_to_ras = diag([1 1 1 1]);
vox_to_ras = vox_to_ras(:);
fwrite(fid,vox_to_ras,'float32');
reserved = zeros(444,1);
fwrite(fid,reserved,'char');
voxel_order = [uint8('RAS')   0];

fwrite(fid,voxel_order,'char');
pad2 = [76     80      83      0];
fwrite(fid,pad2,'char');
image_orientation_patient = [1 0  0   0 1 0];
fwrite(fid,image_orientation_patient,'float32');
pad1 = [0 0];
fwrite(fid,pad1,'char');
invert_xyz = zeros(3,1);
fwrite(fid,invert_xyz,'uchar');
swap_xyz = zeros(3,1);
fwrite(fid,swap_xyz,'uchar');
n_count = length(LineInd);
fwrite(fid,n_count,'int32');
version = 1;
fwrite(fid,version,'int32');
hdr_size = 1000;
fwrite(fid,hdr_size,'int32');

count = 0;
PtVec = PtVec';
PtVec = PtVec(:);
for i = 1:length(LineInd)
    if i==1
        num_pt = LineInd(i);
    else
        num_pt = LineInd(i)-LineInd(i-1);
    end;
    fwrite(fid,num_pt,'int32');
    fwrite(fid,PtVec(3*count+1:3*LineInd(i)),'float32');
    count = LineInd(i);
end;
fclose(fid);
