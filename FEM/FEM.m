function FEM()
    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    z_max = pi;
    r_max = 2;
    h_density = 0.15
    pfix=[0,0;z_max,cos(z_max)]
    [V,F]=distmesh2d(inline('dfct(rz,@(rz) sin(rz(:,1)))','rz'),@huniform,h_density,[0,0;z_max,r_max],[]);
    V;
    F;
end