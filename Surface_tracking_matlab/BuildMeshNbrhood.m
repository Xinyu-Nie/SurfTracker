function [Vert2FaceMap,TrgNbrMap] = BuildMeshNbrhood(coord,trg)
%[Vert2FaceMap,TrgNbrMap] = BuildMeshNbrhood(coord,trg)
%
% Given a triangular mesh, find the neiboring triangles for each vertex and
% triangle. This enables us to access neighboring elements in the mesh.
%
% Inputs:
% (coord,trg): the mesh
%
% Outputs: 
% Vert2FaceMap: index to neighboring triangles of each vertex. The first
% element of each row stores the number of neighboring triangles and the
% rest of the row stores the triangle indexes.
%
% TrgNbrMap: for each triangle, store its three triangle neighbors. 
%
% Yonggang Shi
% Sep 29, 2011.


Vert2FaceMap = zeros(size(coord,1),20);
Vert2FaceMap(:,1) = 1;
for i=1:size(trg,1)
   j =1;
   Vert2FaceMap(trg(i,j),1) = Vert2FaceMap(trg(i,j),1)+1;    
   Vert2FaceMap(trg(i,j),Vert2FaceMap(trg(i,j),1)) = i;
   j =2;
   Vert2FaceMap(trg(i,j),1) = Vert2FaceMap(trg(i,j),1)+1;    
   Vert2FaceMap(trg(i,j),Vert2FaceMap(trg(i,j),1)) = i;
   j =3;
   Vert2FaceMap(trg(i,j),1) = Vert2FaceMap(trg(i,j),1)+1;    
   Vert2FaceMap(trg(i,j),Vert2FaceMap(trg(i,j),1)) = i;   
end;

TrgNbrMap = zeros(size(trg,1),3);

for i=1:size(trg,1)
    %faces sharing vert1 and 2
    for j=2:Vert2FaceMap(trg(i,1),1)
        for k=2:Vert2FaceMap(trg(i,2),1)
            if Vert2FaceMap(trg(i,1),j)~=i && Vert2FaceMap(trg(i,1),j)== Vert2FaceMap(trg(i,2),k)                 
                TrgNbrMap(i,3) = Vert2FaceMap(trg(i,1),j);
                break;
            end;
        end;
        if TrgNbrMap(i,3)>0
            break;
        end;
    end;
    
    
    %faces sharing vert1 and 3
     for j=2:Vert2FaceMap(trg(i,1),1)
        for k=2:Vert2FaceMap(trg(i,3),1)
            if Vert2FaceMap(trg(i,1),j)~=i && Vert2FaceMap(trg(i,1),j)== Vert2FaceMap(trg(i,3),k)                 
                TrgNbrMap(i,2) = Vert2FaceMap(trg(i,1),j);
                break;
            end;
        end;
        if TrgNbrMap(i,2)>0
            break;
        end;
    end;
    %faces sharing vert2 and 3
     for j=2:Vert2FaceMap(trg(i,2),1)
        for k=2:Vert2FaceMap(trg(i,3),1)
            if Vert2FaceMap(trg(i,2),j)~=i && Vert2FaceMap(trg(i,2),j)== Vert2FaceMap(trg(i,3),k)                 
                TrgNbrMap(i,1) = Vert2FaceMap(trg(i,2),j);
                break;
            end;
        end;
        if TrgNbrMap(i,1)>0
            break;
        end;
    end;
end;

