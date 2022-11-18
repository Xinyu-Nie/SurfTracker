function TrackVis2MRTrixTracks(TrackVisFile, NII_TemplateFile, OutputMRTrixTrackFile)
%
%TrackVis2MRTrixTracks(TrackVisFile, NII_TemplateFile, OutputMRTrixTrackFile)
%

[LineInd,PtVec,dim,voxel_size] = ReadTrackVisFile(TrackVisFile);

% tracks = read_mrtrix_tracks(MRTrixTemplateFile);
nii = load_untouch_nii(NII_TemplateFile);

R = [nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z];
sz = nii.hdr.dime.dim(2:4);
s = nii.hdr.dime.pixdim(2:4);

T = R(:,4);
R = R(:,1:3);
R = R*diag(1./s); 

%put back to physical space
PtVec = (R*PtVec')';
for i = 1:3
    PtVec(:,i) = PtVec(:,i) + T(i);
end

%assign to the tracks.data
data = cell(length(LineInd),1);

for i = 1:length(LineInd)
    if i==1
        data{i} = PtVec(1:LineInd(i),:);
    else
        data{i} = PtVec(LineInd(i-1)+1:LineInd(i),:);
    end;
end;

tracks.data = data;
tracks.count = length(LineInd);
tracks.total_count = LineInd(end);

write_mrtrix_tracks(tracks,OutputMRTrixTrackFile);
