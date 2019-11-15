function ComputeEigsTE()
    clf
    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    %create basic rect mesh
    Nz=3;
    Nr=3;
    [V, F] = meshRectangle([0,0],[1,1], Nz, Nr);
    V = V(:,1:2); %vertices blijkbaar 3D punten...
    viewMesh(V, F);
    %fancy algemene boundary mesh
    %z_max = pi;
    %r_max = 1;
    %vertices_length = 0.15;
    %coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
    %[V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    
    for n = 1:length(F)
        %make triangles positively oriented (in de praktijk maken distmesh en
        %meshRect al pos georienteerd dus vrij nutteloos allemaal)
        %if the 3rd point is to the left of the vector point 1-->2, flip the 3rd point with the 2nd
        isLeft = sign((V(F(n,3),1) - V(F(n,1),1)) * (V(F(n,2),2) - V(F(n,1),2)) - (V(F(n,3),2) - V(F(n,1),2)) * (V(F(n,2),1) - V(F(n,1),1)));
        if isLeft > 0
            isLeft
            F(n,[2, 3]) = F(n,[3, 2]);
        elseif isLeft == 0
            disp('error with changing triangle node indices')
        end
        
        z = [V(F(n,1),1), V(F(n,2),1), V(F(n,3),1)];
        r = [V(F(n,1),2), V(F(n,2),2), V(F(n,3),2)];
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
    end
end 