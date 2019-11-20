function ComputeEigsTE()
    clf
    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    %create basic rect mesh
    %%{
    z_max = 1;
    r_max = 1;
    Nz=3;
    Nr=3;
    [V, F] = meshRectangle([0,0],[z_max,r_max], Nz, Nr);
    V = V(:,1:2); %vertices blijkbaar 3D punten...
    %%}
    
    %fancy general boundary mesh
    %{
    z_max = pi;
    r_max = 1;
    vertices_length = 0.15;
    coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
    [V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    %}
    
    % THE SLOW WAY
    %{ 
    %tic
    %pre-allocate to increase performance
    w_nodes = cell(length(F),3); %matlab can't have function handles in array so use cells
    w_edges_z = cell(length(F),3); 
    w_edges_r = cell(length(F),3); 
    w_edges_z_total = cell(length(F),1);
    w_edges_r_total = cell(length(F),1);
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
        %pre-allocate to increase performance
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
            w_edges_z{n,i+1}=@(z,r) 1/(4*A^2)*(a(mod(i+1,3)+1)*b(mod(i+2,3)+1)-a(mod(i+2,3)+1)*b(mod(i+1,3)+1)+(c(mod(i+1,3)+1)*b(mod(i+2,3)+1)-c(mod(i+2,3)+1)*b(mod(i+1,3)+1))*r);
            w_edges_r{n,i+1}=@(z,r) 1/(4*A^2)*(a(mod(i+1,3)+1)*c(mod(i+2,3)+1)-a(mod(i+2,3)+1)*c(mod(i+1,3)+1)+(b(mod(i+1,3)+1)*c(mod(i+2,3)+1)-b(mod(i+2,3)+1)*c(mod(i+1,3)+1))*z);
        end
        %Totale w_edges. Geef in viewMeshandasisFcts altijd point=1 mee
        w_edges_z_total{n,1}=@(z,r) w_edges_z{n,1}(z,r)+w_edges_z{n,2}(z,r)+w_edges_z{n,3}(z,r);
        w_edges_r_total{n,1}=@(z,r) w_edges_r{n,1}(z,r)+w_edges_r{n,2}(z,r)+w_edges_r{n,3}(z,r);
    end
    %toc
    %}
    
    % THE FAST WAY 
    %%{
    %tic
    F=F.';
    V=V.';
    z=[V(1,F(1,:));V(1,F(2,:));V(1,F(3,:))];
    r=[V(2,F(1,:));V(2,F(2,:));V(2,F(3,:))];
    %pre-allocate to increase performance
    a = zeros(3,length(F)); b = zeros(3,length(F)); c = zeros(3,length(F)); A=zeros(1,length(F));
    w_nodes = cell(length(F),3); %matlab can't have function handles in array so use cells
    w_edges_z = cell(length(F),3); w_edges_r = cell(3,length(F)); 
    w_edges_z_total = cell(length(F),1); w_edges_r_total = cell(length(F),1);

    a([1,2,3],:)=z(mod([1,2,3],3)+1,:).*r(mod([1,2,3]+1,3)+1,:)-z(mod([1,2,3]+1,3)+1,:).*r(mod([1,2,3],3)+1,:);
    b([1,2,3],:)=r(mod([1,2,3],3)+1,:)-r(mod([1,2,3]-2,3)+1,:);
    c([1,2,3],:)=z(mod([1,2,3]-2,3)+1,:)-z(mod([1,2,3],3)+1,:);
    A(1,:)=(b(2,:).*c(3,:)-b(3,:).*c(2,:))/2; %idem for all i
    
    for i = 0:2
        for n=1:length(F)
            w_nodes{n,i+1}=@(z,r) (a(i+1,n)+b(i+1,n)*z+c(i+1,n)*r)./(2*A(1,n));
            w_edges_z{n,i+1}=@(z,r) (a(mod(i+1,3)+1,n).*b(mod(i+2,3)+1,n)-a(mod(i+2,3)+1,n).*b(mod(i+1,3)+1,n)+(c(mod(i+1,3)+1,n).*b(mod(i+2,3)+1,n)-c(mod(i+2,3)+1,n).*b(mod(i+1,3)+1,n))*r)./(4*A(1,n).^2);
            w_edges_r{n,i+1}=@(z,r) (a(mod(i+1,3)+1,n).*c(mod(i+2,3)+1,n)-a(mod(i+2,3)+1,n).*c(mod(i+1,3)+1,n)+(b(mod(i+1,3)+1,n).*c(mod(i+2,3)+1,n)-b(mod(i+2,3)+1,n).*c(mod(i+1,3)+1,n))*z)./(4*A(1,n).^2);
        end
    end
    w_edges_z_total{n,1}=@(z,r) w_edges_z{n,1}(z,r)+w_edges_z{n,2}(z,r)+w_edges_z{n,3}(z,r);
    w_edges_r_total{n,1}=@(z,r) w_edges_r{n,1}(z,r)+w_edges_r{n,2}(z,r)+w_edges_r{n,3}(z,r);
    %toc
    %%}
    
    Ig=F([1,2,3,1,2,3,1,2,3],:);
    Jg=F([1,1,1,2,2,2,3,3,3],:);
    Kg=zeros(9,length(F));
    Kg(1,:)=10000;
    Kg(2,:)=100;
    Kg(3,:)=1;
    Kg(5,:)=10000;
    Kg(6,:)=100;
    Kg(9,:)=10000;
    Kg([4,7,8],:)=Kg([2,3,6],:);
    M=full(sparse(Ig,Jg,Kg,length(V),length(V)))
    V=V.';
    F=F.';
    
    %basis rectangle mesh
    %viewMeshandBasisFcts(V, F, w_nodes, w_edges_z, w_edges_r, 1, 1, 1/(6*(Nz+Nr)/2),z_max,r_max)
    %fancy general mesh
    %viewMeshandBasisFcts(V, F, w_nodes, w_edges_z, w_edges_r, 14, 2, vertices_length/20,z_max,r_max)
end 