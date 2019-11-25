function ComputeEigsTE()
    clf
    close all
    addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    %create basic rect mesh
    %%{
    z_max = 1; r_max = 1;
    Nz=15; Nr=Nz;
    [V, F] = meshRectangle([0,0],[z_max,r_max], Nz, Nr);
    V = V(:,1:2); %vertices blijkbaar 3D punten...
    %%}
    
    %fancy general boundary mesh
    %{
    z_max = pi; r_max = 1;
    vertices_length = 0.15;
    coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
    [V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    %}
    
    % THE SLOW WAY
    %{ 
    tic
    %pre-allocate to increase performance
    w_nodes = cell(length(F),3); %matlab can't have function handles in array so use cells
    w_edges_z = cell(length(F),3); 
    w_edges_r = cell(length(F),3); 
    w_edges_z_total = cell(length(F),1);
    w_edges_r_total = cell(length(F),1);
    for n = 1:length(F)
        %make triangles positively oriented (in practice, the mesh algorithm already does this)
        %{
        %if the 3rd point is to the left of the vector point 1-->2, flip the 3rd point with the 2nd
        isLeft = sign((V(F(n,3),1) - V(F(n,1),1)) * (V(F(n,2),2) - V(F(n,1),2)) - (V(F(n,3),2) - V(F(n,1),2)) * (V(F(n,2),1) - V(F(n,1),1)));
        if isLeft > 0
            disp('had to adjust orientation of a triangle!')
            F(n,[2, 3]) = F(n,[3, 2]);
        elseif isLeft == 0
            disp('error with changing triangle node indices')
        end
        %}
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
    toc
    %}
    
    % THE FAST WAY 
    %%{
    tic
    F=F.';
    V=V.';
    z=[V(1,F(1,:));V(1,F(2,:));V(1,F(3,:))];
    r=[V(2,F(1,:));V(2,F(2,:));V(2,F(3,:))];
    %pre-allocate to increase performance
    a = zeros(3,length(F)); b = zeros(3,length(F)); c = zeros(3,length(F)); A=zeros(1,length(F));
    w_nodes = cell(length(F),3); %matlab can't have function handles in arrays so use cells
    w_edges_z = cell(length(F),3); w_edges_r = cell(3,length(F)); 
    w_edges_z_total = cell(length(F),1); w_edges_r_total = cell(length(F),1);
    i=[1,2,3]; s=@(i) mod(i-1,3)+1;
    
    a(i,:)=z(s(i+1),:).*r(s(i+2),:)-z(s(i+2),:).*r(s(i+1),:);
    b(i,:)=r(s(i+1),:)-r(s(i-1),:);
    c(i,:)=z(s(i-1),:)-z(s(i+1),:);
    A(i,:)=(b(s(i+1),:).*c(s(i+2),:)-b(s(i+2),:).*c(s(i+1),:))/2; %idem for all i
    
    for i = 1:3
        for n=1:length(F)
            w_nodes{n,i}=@(z,r) (a(s(i),n)+b(s(i),n)*z+c(s(i),n)*r)./(2*A(i,n));
            w_edges_z{n,i}=@(z,r) (a(s(i+1),n).*b(s(i+2),n)-a(s(i+2),n).*b(s(i+1),n)+(c(s(i+1),n).*b(s(i+2),n)-c(s(i+2),n).*b(s(i+1),n))*r)./(4*A(i,n).^2);
            w_edges_r{n,i}=@(z,r) (a(s(i+1),n).*c(s(i+2),n)-a(s(i+2),n).*c(s(i+1),n)+(b(s(i+1),n).*c(s(i+2),n)-b(s(i+2),n).*c(s(i+1),n))*z)./(4*A(i,n).^2);
        end
    end

    i=[1,2,3];
    w_edges_z_total{n,1}=@(z,r) w_edges_z{n,1}(z,r)+w_edges_z{n,2}(z,r)+w_edges_z{n,3}(z,r);
    w_edges_r_total{n,1}=@(z,r) w_edges_r{n,1}(z,r)+w_edges_r{n,2}(z,r)+w_edges_r{n,3}(z,r);
    toc
    %%}
    
    I=F([1,2,3,1,2,3,1,2,3],:); J=F([1,1,1,2,2,2,3,3,3],:);
    
    m_edge=zeros(9,length(F));
    m_edge([1:9],:)=ones(9,1).*(r(1,:)+r(2,:)+r(3,:))./(2*A(1,:)); 
    M_edge=full(sparse(I,J,m_edge,length(V),length(V)));
    figure('Name','M_edge')
    spy(M_edge)
    
    g_edge=zeros(9,length(F));
    g_edge([1,5,9],:)=((b(s(i+2),:).^2+c(s(i+2),:).^2).*(2*r(s(i),:)+6*r(s(i+1),:)+2*r(s(i+2),:))-2*(b(s(i+1),:).*b(s(i+2),:)+c(s(i+1),:).*c(s(i+2),:)).*(r(s(i),:)+2*r(s(i+1),:)+2*r(s(i+2),:))+(b(s(i+1),:).^2+c(s(i+1),:).^2).*(2*r(s(i),:)+2*r(s(i+1),:)+6*r(s(i+2),:)))./(240*A(i,:));
    g_edge([2,3,6],:)=((b(s(i+2),:).*b(s(i),:)+c(s(i+2),:).*c(s(i),:)).*(r(s(i),:)+2*r(s(i+1),:)+2*r(s(i+2),:))-(b(s(i+2),:).^2+c(s(i+2),:).^2).*(2*r(s(i),:)+2*r(s(i+1),:)+r(s(i+2),:))-(b(s(i+1),:).*b(s(i),:)+c(s(i+1),:).*c(s(i),:)).*(2*r(s(i),:)+2*r(s(i+1),:)+6*r(s(i+2),:))+(b(s(i+1),:).*b(s(i+2),:)+c(s(i+1),:).*c(s(i+2),:)).*(2*r(s(i),:)+r(s(i+1),:)+2*r(s(i+2),:)))./(240*A(i,:))./(240*A(i,:));
    g_edge([4,7,8],:)=g_edge([2,3,6],:); %interaction integrals are symmetric in a triangle's edge indices
    G_edge=sparse(I,J,g_edge,length(V),length(V));
    figure('Name','G_edge')
    spy(G_edge)
    
    g_node=zeros(9,length(F));
    g_node([1,5,9],:)=(3*r(s(i),:)+r(s(i+1),:)+r(s(i+2),:))./(30*A(1,:));
    g_node([2,3,6],:)=(2*r(s(i),:)+2*r(s(i+1),:)+r(s(i+2),:))./(60*A(1,:));
    g_node([4,7,8],:)=g_node([2,3,6],:); %interaction integrals are symmetric in a triangle's edge indices
    G_node=sparse(I,J,g_node,length(V),length(V));
    figure('Name','G_node')
    spy(G_node)
    
    [E,D]=eig(full(G_edge),full(M_edge));
    figure()
    imagesc(reshape(E(1,:),[Nz,Nz]))
    colorbar
    figure()
    imagesc(reshape(E(2,:),[Nz,Nz]))
    colorbar
    figure()
    imagesc(reshape(E(3,:),[Nz,Nz]))
    colorbar
    
    V=V.';
    F=F.';
    
    %basic rectangle mesh
    if (Nz+Nr)/2 < 10
        figure('Name','Mesh and basis functions')
        viewMeshandBasisFcts(V, F, w_nodes, w_edges_z, w_edges_r, 1, 1, 1/(6*(Nz+Nr)/2),z_max,r_max)
    else
        fprintf('are you sure you want to display so many triangles? (%d,%d) \n', Nz, Nr)
    end
    %fancy general mesh
    %{
    if 20/vertices_length < 10
        viewMeshandBasisFcts(V, F, w_nodes, w_edges_z, w_edges_r, 14, 2, vertices_length/20,z_max,r_max)
    else
        fprintf('are you sure you want to display so many triangles? (%d)', 20/vertices_length)
    end
    %}
end 