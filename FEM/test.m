Nz = 5;
Nr = 5;
[Vertices, Faces] = meshRectangle([0,0], [1,1], Nz, Nr);
viewMesh(Vertices(:,1:2), Faces);
