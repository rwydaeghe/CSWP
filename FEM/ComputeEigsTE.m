function ComputeEigsTE()
    clf
    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    %create basic rect mesh
    Nz=3;
    Nr=3;
    [V, F] = meshRectangle([0,0],[1,1], Nz, Nr);
    V = V(:,1:2); %vertices blijkbaar 3D punten...

    %fancy algemene boundary mesh
    %z_max = pi;
    %r_max = 1;
    %vertices_length = 0.05;
    %coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
    %[V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    
    w_nodes = cell(length(F),3); %this special array can contain function handles
    w_edges = cell(length(F),3); 
    w_edges_total = cell(length(F));
    for n = 1:length(F)
        %make triangles positively oriented (in practice, the mesh algorithm already does this)
        %if the 3rd point is to the left of the vector point 1-->2, flip the 3rd point with the 2nd
        isLeft = sign((V(F(n,3),1) - V(F(n,1),1)) * (V(F(n,2),2) - V(F(n,1),2)) - (V(F(n,3),2) - V(F(n,1),2)) * (V(F(n,2),1) - V(F(n,1),1)));
        if isLeft > 0
            disp('had to adjust orientation of a triangle!')
            F(n,[2, 3]) = F(n,[3, 2]);
        elseif isLeft == 0
            disp('error with changing triangle node indices')
        end
        
        z = [V(F(n,1),1); V(F(n,2),1); V(F(n,3),1)];
        r = [V(F(n,1),2); V(F(n,2),2); V(F(n,3),2)];
        a = zeros(3,1);
        b = zeros(3,1);
        c = zeros(3,1);
        for i = 0:2 %moet ik doen want indices in matlab sucken echt
            a(i+1)=z(mod(i+1,3)+1)*r(mod(i+2,3)+1)-z(mod(i+2,3)+1)*r(mod(i+1,3)+1);
            b(i+1)=r(mod(i+1,3)+1)-r(mod(i-1,3)+1);
            c(i+1)=z(mod(i-1,3)+1)-z(mod(i+1,3)+1);
        end
        i=0; %A is toch onafh van i
        A=(b(mod(i+1,3)+1)*c(mod(i+2,3)+1)-b(mod(i+2,3)+1)*c(mod(i+1,3)+1))/2;
        for i = 0:2 
            w_nodes{n,i+1}=@(z,r) 1/(2*A)*(a(i+1)+b(i+1)*z+c(i+1)*r);
            w_edges{n,i+1}=@(z,r) 1/(4*A^2)*[a(mod(i+1,3)+1)*b(mod(i+2,3)+1)-a(mod(i+2,3)+1)*b(mod(i+1,3)+1)+(c(mod(i+1,3)+1)*b(mod(i+2,3)+1)-c(mod(i+2,3)+1)*b(mod(i+1,3)+1))*r;
                                             a(mod(i+1,3)+1)*c(mod(i+2,3)+1)-a(mod(i+2,3)+1)*c(mod(i+1,3)+1)+(b(mod(i+1,3)+1)*c(mod(i+2,3)+1)-b(mod(i+2,3)+1)*c(mod(i+1,3)+1))*z];
        end
        %To see total w_edges in a triangle, and edit arguments in viewMesh
        %w_edges_total{n}=@(z,r) w_edges{n,1}(z,r)+w_edges{n,2}(z,r)+w_edges{n,3}(z,r);
    end
    viewMeshandBasisFcts(V, F, w_nodes, w_edges, 1, 1, 1/(6*(Nz+Nr)/2))
end 