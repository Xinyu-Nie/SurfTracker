function [pt_coordinates,pt_normals,triangles,pt_colors] = ReadObjShape(fname)
% [pt_coordinates,pt_normals,triangles] = ReadObjShape(fname)
% This function reads a shape of OBJ format. 
%
%Yonggang Shi
%May 2015

fid = fopen(fname,'r');
shape_type = fscanf(fid,'%s',1);

if shape_type == 'P' %Polygon
    a = fscanf(fid,'%f',6);
    num_pts = a(6);
    pt_coordinates = fscanf(fid,'%f',num_pts*3);
    pt_coordinates = reshape(pt_coordinates,3,num_pts);
    pt_coordinates = pt_coordinates';
    
    pt_normals = fscanf(fid,'%f',num_pts*3);
    pt_normals = reshape(pt_normals,3,num_pts);
    pt_normals = pt_normals';

    num_triangles = fscanf(fid,'%f',1);
    color_type = fscanf(fid,'%f',1);
    if color_type==2  
        pt_colors = fscanf(fid,'%f',num_pts*4);
        pt_colors = reshape(pt_colors,4,num_pts);
        pt_colors = pt_colors';
    elseif color_type==0 
        pt_colors = fscanf(fid,'%f',4);%num_pts*4);
    end;

    tmp = fscanf(fid,'%f',num_triangles);
    triangles = fscanf(fid,'%f',num_triangles*3);  
    triangles = triangles +1; %adjust for matlab index starting from 1 instead of 0
    triangles = reshape(triangles,3, num_triangles);
    triangles = triangles';
        
end;

if shape_type == 'L' %line
    a = fscanf(fid,'%f',2);
    num_pts = a(2);
    pt_coordinates = fscanf(fid,'%f',num_pts*3);
    pt_coordinates = reshape(pt_coordinates,3,num_pts);
    pt_coordinates = pt_coordinates';

end;

fclose(fid);
