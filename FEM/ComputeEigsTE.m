function ComputeEigsTE(meshSize)
    global N
    N=meshSize; 
    %N=5; 
    clf; close all; addpath('./DistMesh'); %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    
    %% Create and analyze mesh
    %create V and F
    meshType='cylinder';
    z_max = 1; r_max = 1;    
    if meshType=='cylinder'
        Nz=N; Nr=Nz;
        global V; global F
        [V, F] = meshRectangle([0,0],[z_max,r_max], Nz, Nr); V = V(:,1:2); %vertices blijkbaar 3D punten...
    else
        vertices_length = 0.4;
        coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
        [V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    end
    V=V.'; F=F.';
    F(:,1:2:end)=F([1,3,2],1:2:end); %tang cont
    
    %Create E
    i_index=[1,2,3,1,2,3,1,2,3];j_index=[1,1,1,2,2,2,3,3,3]; I=F(i_index,:); J=F(j_index,:);
    e=zeros(9,size(F,2)); e([2,3,6],:)=1;
    ME1=sparse(I,J,e,length(V),length(V)); 
    E_boundary=find2D(triu(ME1+ME1.')==1); E_internal=find2D(triu(ME1+ME1.')==2); global E; E=find2D(triu(ME1+ME1.'));
    
    %Split E in internal and boundary edges
    if meshType == 'cylinder'
        %Maar O(N^2) want slechts boundary. 3rde grootste bottleneck na Fedge en eigs
        E_vec=round((V(:,imag(E_boundary))-V(:,real(E_boundary)))*(N-1)); %normalized to be logical array. round for floating point error
        E_hori=E_boundary(E_vec(1,:)==1); E_vert=E_boundary(E_vec(2,:)==1);
        E_axis=E_hori(V(2,real(E_hori))==0); E_mantle=E_hori(V(2,real(E_hori))==1); E_leftBound=E_vert(V(1,real(E_vert))==0); E_rightBound=E_vert(V(1,real(E_vert))==1);
        E_boundary=[E_leftBound,E_rightBound,E_axis,E_mantle];
        E_boundary_edge=reshape(full(sum(sparse(bsxfun(@eq,E_boundary(:),E.').*[1:length(E)]),2)),[],4);
    end
    
    %RVWn?
    % neumann is altijd goed (kijk heel goed naar interpretatie d/dz e_z
    % met alle drie w_edges en dan op figuur kijken naar pijltjes steeds
    % zelfde z comp, of bewijs wiskundig). Axis is problematisch want daar moet
    % e_r=0 in TM op edges en het beste dat ik kan doen e_z=0 door
    % horizontale edge =0 te doen. Anders kun je zeggen dat het per default al ok is?
    % Dirichlet mantel is te doen door horizontale weg te doen en dan hebt ge e_z=0 zou gewild
    
    %Create Fedge (F as function of E)
    e=e.*[1:size(F,2)]; %introduces a face label when edges are counted
    global Fedge; Fedge=zeros(3,size(F,2)); Coordn=[I(:),J(:),e(:)]; %labels in third dimension
    ME2=ndSparse.build(Coordn(e(:)~=0,:),1,[length(V),length(V),size(F,2)]);
    Fedge=find3D(findPattern(find3D(ME2+permute(ME2,[2,1,3]),'3d'), E.'),'2d'); %Dit algoritme behoudt de alternerende pos of neg orientatie
    Fedge([1,2,3],2:2:end)=Fedge([3,1,2],2:2:end); %reshape(E(Fedge(:)),3,[])==F(2,:)+j*F(3,:) %(onschuldige) even permutatie zodat Fedge(1,:)<->F(1,:) en bewijs hiervoor  
    
    %% Additional variables. 
    %Pre-allocate to increase performance
    z=[V(1,F(1,:));V(1,F(2,:));V(1,F(3,:))]; r=[V(2,F(1,:));V(2,F(2,:));V(2,F(3,:))];
    a = zeros(3,size(F,2)); b = zeros(3,size(F,2)); c = zeros(3,size(F,2)); A=zeros(1,size(F,2));
    global w_nodes; global w_edges_z; global w_edges_r; global w_edges_z_total; global w_edges_r_total;
    w_nodes = cell(size(F,2),3); %matlab can't have function handles in arrays so use cells
    w_edges_z = cell(size(F,2),3); w_edges_r = cell(size(F,2),3); 
    w_edges_z_total = cell(size(F,2),1); w_edges_r_total = cell(size(F,2),1);
    
    i=[1,2,3]; s=@(i) mod(i-1,3)+1; %cyclical indexing
    a(i,:)=z(s(i+1),:).*r(s(i+2),:)-z(s(i+2),:).*r(s(i+1),:);
    b(i,:)=r(s(i+1),:)-r(s(i-1),:);
    c(i,:)=z(s(i-1),:)-z(s(i+1),:);
    A(i,:)=(b(s(i+1),:).*c(s(i+2),:)-b(s(i+2),:).*c(s(i+1),:))/2; %idem for all i
    %a(:,2:2:end)=-a([1,3,2],2:2:end); b(:,2:2:end)=-b([1,3,2],2:2:end); c(:,2:2:end)=-c([1,3,2],2:2:end);
    
    %%{
    for i = 1:3
        for n=1:size(F,2)
            %idee vectorieel en function handle opchrijven als str. Die dan
            %omzetten naar code wanneer nodig
            %Let op ik ben NIET zeker of ik hier abs(A) moet doen...
            %Nevermind, je haalt dat gewoon uit moet = 1 op het punt
            w_nodes{n,i}=@(z,r) (a(s(i),n)+b(s(i),n)*z+c(s(i),n)*r)./(2*A(i,n));
            w_edges_z{n,i}=@(z,r) (a(s(i+1),n).*b(s(i+2),n)-a(s(i+2),n).*b(s(i+1),n)+(c(s(i+1),n).*b(s(i+2),n)-c(s(i+2),n).*b(s(i+1),n))*r)./(4*A(i,n).^2);
            w_edges_r{n,i}=@(z,r) (a(s(i+1),n).*c(s(i+2),n)-a(s(i+2),n).*c(s(i+1),n)+(b(s(i+1),n).*c(s(i+2),n)-b(s(i+2),n).*c(s(i+1),n))*z)./(4*A(i,n).^2);
            w_edges_z_total{n,1}=@(z,r) w_edges_z{n,1}(z,r)+w_edges_z{n,2}(z,r)+w_edges_z{n,3}(z,r);
            w_edges_r_total{n,1}=@(z,r) w_edges_r{n,1}(z,r)+w_edges_r{n,2}(z,r)+w_edges_r{n,3}(z,r);    
        end
    end
    %}
    i=[1,2,3];
    
    %% Build matrices and find eigenvalues
    I_edge=Fedge(i_index,:); J_edge=Fedge(j_index,:);
    
    m_edge=zeros(9,size(F,2));
    m_edge(1:9,:)=ones(9,1).*(r(1,:)+r(2,:)+r(3,:))./(2*abs(A(1,:)));%let op eigenlijk is A in absolute waarde. Je kunt dat ook nagaan uit symmetrieredenen LL
    M_edge=sparse(I_edge,J_edge,m_edge,length(E),length(E));

    g_edge=zeros(9,size(F,2));
    g_edge([1,5,9],:)=((b(s(i+2),:).^2+c(s(i+2),:).^2).*(2*r(s(i),:)+6*r(s(i+1),:)+2*r(s(i+2),:))-2*(b(s(i+1),:).*b(s(i+2),:)+c(s(i+1),:).*c(s(i+2),:)).*(r(s(i),:)+2*r(s(i+1),:)+2*r(s(i+2),:))+(b(s(i+1),:).^2+c(s(i+1),:).^2).*(2*r(s(i),:)+2*r(s(i+1),:)+6*r(s(i+2),:)))./(240*abs(A(i,:)));
    g_edge([2,6,3],:)=((b(s(i+2),:).*b(s(i),:)+c(s(i+2),:).*c(s(i),:)).*(r(s(i),:)+2*r(s(i+1),:)+2*r(s(i+2),:))-(b(s(i+2),:).^2+c(s(i+2),:).^2).*(2*r(s(i),:)+2*r(s(i+1),:)+r(s(i+2),:))-(b(s(i+1),:).*b(s(i),:)+c(s(i+1),:).*c(s(i),:)).*(2*r(s(i),:)+2*r(s(i+1),:)+6*r(s(i+2),:))+(b(s(i+1),:).*b(s(i+2),:)+c(s(i+1),:).*c(s(i+2),:)).*(2*r(s(i),:)+r(s(i+1),:)+2*r(s(i+2),:)))./(240*abs(A(i,:)));
    %equiv tangential cont: switch in 2<->3 and 5<->9 for tangential continuity every other triangle
    %g_edge([2,3,5,9],2:2:end)=g_edge([3,2,9,5],2:2:end);
    g_edge([4,7,8],:)=g_edge([2,3,6],:); %interaction integrals are symmetric in a triangle's edge indices
    G_edge=sparse(I_edge,J_edge,g_edge,length(E),length(E));
    
    g_node=zeros(9,size(F,2));
    g_node([1,5,9],:)=(3*r(s(i),:)+r(s(i+1),:)+r(s(i+2),:))./(30*abs(A(i,:)));
    g_node([2,6,3],:)=(2*r(s(i),:)+2*r(s(i+1),:)+r(s(i+2),:))./(60*abs(A(i,:)));
    g_node([4,7,8],:)=g_node([2,3,6],:); %interaction integrals are symmetric in a triangle's edge indices
    G_node=sparse(I,J,g_node,length(V),length(V));
    
    m_node=zeros(9,size(F,2));
    m_node(1:9,:)=c(i_index,:)/6+c(j_index,:)/6+(c(i_index,:).*c(j_index,:)+b(i_index,:).*b(j_index,:))./(2*abs(A(i_index,:))).*(r(1,:)+r(2,:)+r(3,:))/6;
    M_node=sparse(I,J,m_node,length(V),length(V));
        
    %RVWn
    %%{
    sz1=size(G_edge,1); sz2=size(E_boundary_edge(:,4),1);
    G_edge(E_boundary_edge(:,4),:)=zeros(sz2,sz1);
    G_edge(:,E_boundary_edge(:,4))=zeros(sz1,sz2);
    M_edge(E_boundary_edge(:,4),:)=zeros(sz2,sz1);
    M_edge(:,E_boundary_edge(:,4))=zeros(sz1,sz2);
    for n=1:size(E_boundary_edge(:,4),1)
        M_edge(E_boundary_edge(n,4),E_boundary_edge(n,4))=1;
        %G_edge(E_boundary_edge(n,4),E_boundary_edge(n,4))=1;
    end
    
    sz1=size(G_edge,1); sz2=size(E_boundary_edge(:,3),1);
    G_edge(E_boundary_edge(:,3),:)=zeros(sz2,sz1);
    G_edge(:,E_boundary_edge(:,3))=zeros(sz1,sz2);
    M_edge(E_boundary_edge(:,3),:)=zeros(sz2,sz1);
    M_edge(:,E_boundary_edge(:,3))=zeros(sz1,sz2);
    for n=1:size(E_boundary_edge(:,3),1)
        M_edge(E_boundary_edge(n,3),E_boundary_edge(n,3))=1;
        %G_edge(E_boundary_edge(n,3),E_boundary_edge(n,3))=1;
    end
    %}

    %length(E)-size(F,2)+1 %eig gwn #nodes volgens euler
    %sum(svd(full(M_edge))<0.01)
    %[Eedge,D]=eig(full(M_edge),500*full(G_edge)); %blijkaar zou dit een juister resultaat geven dan eigs
    myeigs=eig(full(M_edge),full(G_edge));
    %sum(myeigs<1)
    %plot([1:1:size(myeigs)],myeigs)
    global Eedge; [Eedge,D]=eig(full(M_edge),full(G_edge)); %also not quite O(n)... %(5*pi)^2=250 en we willen geen complexe getallen
    dnnz=sort(D(D>1)); 
    %plot([1:1:size(sort(myeigs))],sort(myeigs))
    %dnnz(1:5)
    %sparse(D./max(D(:)))
    
    %%{
    if meshType == 'cylinder'
        if (Nz+Nr)/2 < 11
            visualizeElectricField(1/(6*(Nz+Nr)/2)/3,z_max,r_max)
            viewMeshandBasisFcts(1, 1, 1/(6*(Nz+Nr)/2),z_max,r_max)
            viewMeshandBasisFcts(2, 1, 1/(6*(Nz+Nr)/2),z_max,r_max)
        else
            fprintf('are you sure you want to display so many triangles? (%d,%d) \n', Nz, Nr)
        end
    else
        if size(F,2) < 21
            viewMeshandBasisFcts(14, 2, vertices_length/20,z_max,r_max)
        else
            fprintf('are you sure you want to display so many triangles? (%d) \n', size(F,2))
        end
    end
    %}
end

function out=trianglesImIn(z,r) 
    global F; global w_nodes; 
    out=[]; 
    for n=1:size(F,2)
        if (w_nodes{n,1}(z,r) <= 1+eps & w_nodes{n,1}(z,r) >= 0-eps) & (w_nodes{n,2}(z,r) <= 1+eps & w_nodes{n,2}(z,r) >= 0-eps) & (w_nodes{n,3}(z,r) <= 1+eps & w_nodes{n,3}(z,r) >= 0-eps)
           out=[out, n];
        end
    end
end

function out=Enodes(Z,R)
    global w_nodes
    out=zeros(size(Z,1),size(R,2));
    for iz=1:size(Z,1)
        for ir=1:size(R,2)
            z=Z(1,iz); r=R(ir,1);
            for tri=trianglesImIn(z,r)
                for node=1:3
                    out(iz,ir)=out(iz,ir)+w_nodes{tri,node}(z,r);
                end
            end
        end
    end    
end

function [Ez,Er]=ElectricField_TM(Z,R)
    global w_edges_z; global w_edges_r; 
    Ez=zeros(size(Z,1),size(R,2)); Er=zeros(size(Z,1),size(R,2));
    for iz=1:size(Z,2)
        for ir=1:size(R,1)
            z=Z(1,iz); r=R(ir,1); 
            for tri=trianglesImIn(z,r) %brute force but visualization is only meant for low N...
                for edge=1:3
                    Ez(iz,ir)=Ez(iz,ir)+w_edges_z{tri,edge}(z,r);
                    Er(iz,ir)=Er(iz,ir)+w_edges_r{tri,edge}(z,r);                    
                end
            end
        end
    end    
end

function out = find3D(X,shape_set) 
    [~,y,z]=size(X);
    flat=mod(find(permute(ndSparse(X),[2,1,3]))-1,y)+1; %let op: doet al een triu (idk hoe?)
    if shape_set == '3d'
        out=reshape(flat(1:2:end)+j*flat(2:2:end),[3,1,z]);
    elseif shape_set == '2d'
        out=reshape(flat,[3,z]);
    end
end

function out = find2D(X) %for elegance in one line and efficient double-element storage
    [x,y]=find(X); out = (x+j*y);
end

function out=findPattern(A,B) %find patternS A in data B and put them in array
    global F; lenF=size(F,2); global N; 
    %bedoeling: O(n) maken van out=bsxfun(@eq,A,B) 
    if N <=7
        out=bsxfun(@eq,A,B); %algoritmes werken eigenlijk enkel voor grote N. Bovendien is dit minder overlay = sneller
    else
        %create O(1) windows where pattern is probable to be found
        axis=[1:length(B)/lenF:length(B)];
        mean2=axis-[1:1:lenF]/lenF*N+ones(1,lenF)*N;
        mean1=mean2-ones(1,lenF)*(3.045*N-3.88); knik=floor(3*N*lenF/length(B)); mean1=[axis(1:knik)/axis(knik)*mean1(knik) mean1(knik+1:end)];
        if N>=50
            halfwindowsize=N/10;
        else
            halfwindowsize=5;
        end
        if N>=41
            %another speed hack for large N
            window1start=mean1; window1start(window1start<1)=1; window1start(window1start>length(B))=length(B); window1start=ceil(window1start);
        else
            window1start=mean1-halfwindowsize; window1start(window1start<1)=1; window1start(window1start>length(B))=length(B); window1start=ceil(window1start);
        end
        window1end=mean1+halfwindowsize; window1end(window1end<1)=1; window1end(window1end>length(B))=length(B); window1end=ceil(window1end);
        window2start=mean2-halfwindowsize; window2start(window2start<1)=1; window2start(window2start>length(B))=length(B); window2start=ceil(window2start);
        window2end=mean2+halfwindowsize; window2end(window2end<1)=1; window2end(window2end>length(B))=length(B); window2end=ceil(window2end);
        coordwindow1=[]; coordwindow2=[]; 
        for u=1:lenF %not too bad since operation itself is O(1)
            newbit1=window1start(u):window1end(u);
            coordwindow1=[coordwindow1 [newbit1;u*ones(size(newbit1))]];
            newbit2=window2start(u):window2end(u);
            coordwindow2=[coordwindow2 [newbit2;u*ones(size(newbit2))]];
        end    
        coordwindows=[coordwindow1 coordwindow2]; sz=size(coordwindows); coordwindows=[coordwindows; ones(1,sz(2))];
        
        %findPattern in windows using sparse matrices
        B=ndSparse.build(coordwindows([3,1,2],:).', B(1,coordwindows(1,:)).');
        out=ndSparse(bsxfun(@eq,A,B));
        %Visualisatie
        %{
        plot(axis,[mean1;window1start;window1end;mean2;window2start;window2end])
        hold on
        out_vis=mod(find(sum(out,1))-1,length(B))+1;    
        scatter([1:length(B)/(length(out_vis)+1):length(B)],out_vis.','.')
        %}
    end
end

function visualizeElectricField(resolution,z_max,r_max)
    global V; global F;
    [Z,R] = meshgrid(0:resolution:z_max,0:resolution:r_max);
    figure('Name','Electric field')
    hold on
    %contourf(Z,R,Enodes(Z,R),1/(5*resolution),'--','ShowText','on');
    V=V.'; F=F.'; patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F,'FaceColor','none','EdgeColor','black','LineWidth',1.5); V=V.'; F=F.';     
    [Ez,Er]=ElectricField_TM(Z,R);
    quiver(Z,R,Ez.',Er.','black'); %let op transpose van matrices omdat meshgrid nogal vreemd is... Kijk naar Z en R

    hold off
end

function viewMeshandBasisFcts(triangle,point,resolution,z_max,r_max)
    global V; global F; global w_nodes; global w_edges_z; global w_edges_r; global w_edges_z_total; global w_edges_r_total;
    V=V.'; %to do: write viewMesh so this isn't necessary(?)
    F=F.';
    [Z,R] = meshgrid(0:resolution:z_max,0:resolution:r_max);
    
    figure('Name','Mesh and basis functions')
    hold on
    %show w_node as filled contourplot
    [c,h] = contourf(Z,R,w_nodes{triangle,point}(Z,R),1/(2*resolution),'--','ShowText','on');
    h.LevelList=round(h.LevelList,2);  %rounds levels to 2nd decimal place
    clabel(c,h);
    %show w_edge as vectorplot
    quiver(Z,R,w_edges_z_total{triangle,point}(Z,R),w_edges_r_total{triangle,point}(Z,R),'black'); 
    %show mesh
    patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F,'FaceColor','none','EdgeColor','black','LineWidth',1.5);
    %highlight selected triangle in blue
    patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F(triangle,:),'FaceColor','none','EdgeColor','blue', 'LineWidth', 1.5);
    %highlight selected node in red
    scatter(V(F(triangle,point),1), V(F(triangle,point),2),'filled','red');        
    %highlight selected edge in red
    p2=mod((point+1)-1,3)+1; p3=mod((point+2)-1,3)+1;
    plot([V(F(triangle,p2),1), V(F(triangle,p3),1)],[V(F(triangle,p2),2), V(F(triangle,p3),2)],'red','LineWidth',1.5);
    xlabel('$z$ [$m$]','interpreter','latex')
    ylabel('$r$ [$m$]','interpreter','latex')
    hold off
    
    V=V.';
    F=F.';
end