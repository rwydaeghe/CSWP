function FEM()
    %Nz=17;
    %Nr=17;
    %[Vertices, Faces] = meshRectangle([0,0],[1,1], Nz, Nr);
    %viewMesh(Vertices, Faces);

    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    fd=@(p) fct(p);
    [p,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
end
function out = fct(p)
    out = p(:,1).^2/2^2+p(:,2).^2/1^2-1;
    %if p(:,1) > 0
    %    out = p(:,1).^2/2^2+p(:,2).^2/1^2-1;
    %else
    %    out = 0;
    %end
end
