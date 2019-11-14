function FEM()
    clf
    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    z_max = pi;
    r_max = 1;
    vertices_length = 0.15;
    coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
    [V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    
    
    
    
    %a=[1,2,3];
    %a(mod(4,3))
end