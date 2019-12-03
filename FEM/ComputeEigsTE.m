function ComputeEigsTE(N)
    %N=180;
    %clf
    %close all
    %addpath('./DistMesh') %Je kunt nu ook scriptjes vinden in die grote folder voor meshes
    %create basic rect mesh
    %%{
    z_max = 1; r_max = 1;
    Nz=N; Nr=Nz;
    [V, F] = meshRectangle([0,0],[z_max,r_max], Nz, Nr);
    V = V(:,1:2); %vertices blijkbaar 3D punten...
    %%}
    
    %fancy general boundary mesh
    %{
    z_max = pi; r_max = 1;
    vertices_length = 0.4;
    coord_fix=[0,0;0,boundary(0);z_max,0;z_max,boundary(z_max)];
    [V,F]=distmesh2d(inline('dfct(rz,@(rz) boundary(rz(:,1)))','rz'),@huniform,vertices_length,[0,0;z_max,r_max],coord_fix);
    %}
    
    V=V.'; F=F.';
    z=[V(1,F(1,:));V(1,F(2,:));V(1,F(3,:))]; r=[V(2,F(1,:));V(2,F(2,:));V(2,F(3,:))];
    %pre-allocate to increase performance
    a = zeros(3,length(F)); b = zeros(3,length(F)); c = zeros(3,length(F)); A=zeros(1,length(F));
    w_nodes = cell(length(F),3); %matlab can't have function handles in arrays so use cells
    w_edges_z = cell(length(F),3); w_edges_r = cell(length(F),3); 
    w_edges_z_total = cell(length(F),1); w_edges_r_total = cell(length(F),1);
    i=[1,2,3]; s=@(i) mod(i-1,3)+1; %cyclical indexing
    
    a(i,:)=z(s(i+1),:).*r(s(i+2),:)-z(s(i+2),:).*r(s(i+1),:);
    b(i,:)=r(s(i+1),:)-r(s(i-1),:);
    c(i,:)=z(s(i-1),:)-z(s(i+1),:);
    A(i,:)=(b(s(i+1),:).*c(s(i+2),:)-b(s(i+2),:).*c(s(i+1),:))/2; %idem for all i
    
    for i = 1:3
        for n=1:length(F)
            %idee vectorieel en function handle opchrijven als str. Die dan
            %omzetten naar code wanneer nodig
            %w_nodes{n,i}=@(z,r) (a(s(i),n)+b(s(i),n)*z+c(s(i),n)*r)./(2*A(i,n));
            %w_edges_z{n,i}=@(z,r) (a(s(i+1),n).*b(s(i+2),n)-a(s(i+2),n).*b(s(i+1),n)+(c(s(i+1),n).*b(s(i+2),n)-c(s(i+2),n).*b(s(i+1),n))*r)./(4*A(i,n).^2);
            %w_edges_r{n,i}=@(z,r) (a(s(i+1),n).*c(s(i+2),n)-a(s(i+2),n).*c(s(i+1),n)+(b(s(i+1),n).*c(s(i+2),n)-b(s(i+2),n).*c(s(i+1),n))*z)./(4*A(i,n).^2);
        end
    end

    i=[1,2,3];
    %w_edges_z_total{n,1}=@(z,r) w_edges_z{n,1}(z,r)+w_edges_z{n,2}(z,r)+w_edges_z{n,3}(z,r);
    %w_edges_r_total{n,1}=@(z,r) w_edges_r{n,1}(z,r)+w_edges_r{n,2}(z,r)+w_edges_r{n,3}(z,r);    
    
    I=F([1,2,3,1,2,3,1,2,3],:); J=F([1,1,1,2,2,2,3,3,3],:);
    
    %Create edge list
    e=zeros(9,length(F)); e([2,3,6],:)=1; e=e.*[1:length(F)];
    ME1=sparse(I,J,e,length(V),length(V)); E=find2D(triu(ME1+ME1.')); %type triu(..)==1,2 to get edge list of resp. once or doubly counted edges
    
    %Write F as function of E      
    Fedge=zeros(3,length(F)); Coordn=[I(:),J(:),e(:)];
    ME2=ndSparse.build(Coordn(e(:)~=0,:),1,[length(V),length(V),length(F)]);
    Fedge=find3D(findPattern(find3D(ME2+permute(ME2,[2,1,3]),'3d'), E.',length(F),N),'2d');
    I=F([1,2,3,1,2,3,1,2,3],:); J=F([1,1,1,2,2,2,3,3,3],:);
    I_edge=Fedge([1,2,3,1,2,3,1,2,3],:); J_edge=Fedge([1,1,1,2,2,2,3,3,3],:);
    
    m_edge=zeros(9,length(F));
    m_edge([1:9],:)=ones(9,1).*(r(1,:)+r(2,:)+r(3,:))./(2*A(1,:));
    %m_edge(:,2:2:end)=-m_edge(:,2:2:end); %tangential continuity
    M_edge=sparse(I_edge,J_edge,m_edge,length(E),length(E));
    %figure('Name','M_edge')
    %spy(M_edge)
    
    g_edge=zeros(9,length(F));
    g_edge([1,5,9],:)=((b(s(i+2),:).^2+c(s(i+2),:).^2).*(2*r(s(i),:)+6*r(s(i+1),:)+2*r(s(i+2),:))-2*(b(s(i+1),:).*b(s(i+2),:)+c(s(i+1),:).*c(s(i+2),:)).*(r(s(i),:)+2*r(s(i+1),:)+2*r(s(i+2),:))+(b(s(i+1),:).^2+c(s(i+1),:).^2).*(2*r(s(i),:)+2*r(s(i+1),:)+6*r(s(i+2),:)))./(240*A(i,:));
    g_edge([2,3,6],:)=((b(s(i+2),:).*b(s(i),:)+c(s(i+2),:).*c(s(i),:)).*(r(s(i),:)+2*r(s(i+1),:)+2*r(s(i+2),:))-(b(s(i+2),:).^2+c(s(i+2),:).^2).*(2*r(s(i),:)+2*r(s(i+1),:)+r(s(i+2),:))-(b(s(i+1),:).*b(s(i),:)+c(s(i+1),:).*c(s(i),:)).*(2*r(s(i),:)+2*r(s(i+1),:)+6*r(s(i+2),:))+(b(s(i+1),:).*b(s(i+2),:)+c(s(i+1),:).*c(s(i+2),:)).*(2*r(s(i),:)+r(s(i+1),:)+2*r(s(i+2),:)))./(240*A(i,:))./(240*A(i,:));
    g_edge([4,7,8],:)=g_edge([2,3,6],:); %interaction integrals are symmetric in a triangle's edge indices
    %g_edge(:,2:2:end)=-g_edge(:,2:2:end); %tangential continuity
    G_edge=sparse(I_edge,J_edge,g_edge,length(E),length(E));
    %figure('Name','G_edge')
    %spy(G_edge)
    
    g_node=zeros(9,length(F));
    g_node([1,5,9],:)=(3*r(s(i),:)+r(s(i+1),:)+r(s(i+2),:))./(30*A(i,:));
    g_node([2,3,6],:)=(2*r(s(i),:)+2*r(s(i+1),:)+r(s(i+2),:))./(60*A(i,:));
    g_node([4,7,8],:)=g_node([2,3,6],:); %interaction integrals are symmetric in a triangle's edge indices
    G_node=sparse(I,J,g_node,length(V),length(V));
    %figure('Name','G_node')
    %spy(G_node)
        
    length(E)-length(F)+1 %eig gwn #nodes volgens euler
    sum(svd(full(M_edge))<0.01)
    [Eedge,D]=eigs(M_edge,G_edge,5); %also not quite O(n)...
    sparse(D)
    %Eedge
    
    V=V.'; %to do: write viewMesh so this isn't necessary(?)
    F=F.';
    
    %basic rectangle mesh
    %%{
    if (Nz+Nr)/2 < 10
        %figure('Name','Mesh and basis functions')
        %viewMeshandBasisFcts(V, F, w_nodes, w_edges_z, w_edges_r, 1, 1, 1/(6*(Nz+Nr)/2),z_max,r_max)
    else
        fprintf('are you sure you want to display so many triangles? (%d,%d) \n', Nz, Nr)
    end
    %%}
    %fancy general mesh
    %{
    if length(F) < 21
        viewMeshandBasisFcts(V, F, w_nodes, w_edges_z, w_edges_r, 14, 2, vertices_length/20,z_max,r_max)
    else
        fprintf('are you sure you want to display so many triangles? (%d) \n', length(F))
    end
    %}
end 

function out = find3D(X,shape_set) 
    [x,y,z]=size(X);
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

function out=findPattern(A,B,lenF,N) %find patternS A in data B and put them in array
    %bedoeling: O(n) maken van out=bsxfun(@eq,A,B) 
    if N <=7
        out=bsxfun(@eq,A,B); %algoritmes werken eigenlijk enkel voor grote N. Bovendien is dit minder overlay = sneller
    else
        %create windows
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
        for u=1:lenF
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