function VertNbrCellArr = MeshNeiboringVertices(coord,trg)
%VertNbrCellArr = MeshNeiboringVertices(coord,trg)
%
%Given a triangular mesh, find the vertices in the one-ring neighborhood of
%each vertex. The output is a cell array. 

[Vert2FaceMap,TrgNbrMap] = BuildMeshNbrhood(coord,trg);
Vert2FaceCellArr = cell(size(coord,1),1);
for i = 1:size(coord,1)
    Vert2FaceCellArr{i} = Vert2FaceMap(i,2:Vert2FaceMap(i,1));
end;

VertNbrCellArr = cell(size(coord,1),1);
for i = 1:size(coord,1)
    a = trg(Vert2FaceCellArr{i},:);
    VertNbrCellArr{i} = setdiff(a(:)',i);
end
