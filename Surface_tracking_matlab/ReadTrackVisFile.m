function [LineInd,PtVec,dim,voxel_size] = ReadTrackVisFile(fname)
%[LineInd,PtVec,dim,voxel_size] = ReadTrackVisFile(fname)
%
% Yonggang Shi
% USC Stevens Neuroimaging and Informatics Institute
% Laboratory of Neuro Imaging (LONI)
% Keck School of Medicine of USC
% University of Southern California


fid = fopen(fname,'r');
id_string = fread(fid,6,'char');
dim = fread(fid,3,'int16');
voxel_size = fread(fid,3,'float32');
origin = fread(fid,3,'float32');
n_s = fread(fid,1,'int16');
s_name = fread(fid,200,'char');
n_p = fread(fid,1,'int16');
p_name = fread(fid,200,'char');
vox_to_ras = fread(fid,16,'float32');
vox_to_ras = reshape(vox_to_ras,[4 4])';
reserved = fread(fid,444,'char');
voxel_order = fread(fid,4,'char');
pad2 = fread(fid,4,'char');
image_orientation_patient = fread(fid,6,'float32');
pad1 = fread(fid,2,'char');
invert_xyz = fread(fid,3,'uchar');
swap_xyz = fread(fid,3,'uchar');
n_count = fread(fid,1,'int32');
version = fread(fid,1,'int32');
hdr_size = fread(fid,1,'int32');

if n_count>0    
    MaxNum = n_count*100;
    PtVec = zeros(MaxNum*(3+n_s),1);
    LineInd = zeros(n_count,1);
    count = 0;
    for i=1:n_count
        num_pt = fread(fid,1,'int32');
        PtVec((3+n_s)*count+1:(3+n_s)*count+(3+n_s)*num_pt) = fread(fid,num_pt*(3+n_s),'float32');
        count = count+num_pt;  
        LineInd(i) = count;    
    end
else
    flag = 1;
    i=0;
    PtVec = [];
    LineInd = [];
    count = 0;
    while flag==1
        [num_pt,num_read] = fread(fid,1,'int32');
        if num_read>0 
            xxx =  fread(fid,num_pt*(3+n_s),'float32');
            if num_pt>100
                PtVec((3+n_s)*count+1:(3+n_s)*count+(3+n_s)*num_pt) = xxx;
                count = count+num_pt;  
                i = i+1;
                LineInd(i) = count;   
            end
        else
            flag = 0;
        end
    end
end
PtVec = PtVec(1:(3+n_s)*count);
PtVec = reshape(PtVec,[(3+n_s) count]);
PtVec = PtVec';
fclose(fid);
