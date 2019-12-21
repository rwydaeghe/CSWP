function ComputeEigsTE(meshSize, mode)
    global N
    N=meshSize; 
    %N=50; 
    %mode=2;    
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
    ECN=sparse(I,J,e,length(V),length(V)); %Edge Count on Node
    E_boundary=find2D(triu(ECN+ECN.')==1); global E; E=find2D(triu(ECN+ECN.')); %note that E's edges don't have correct "directions Re->Im" but not very important
    
    %Create Fedge (F as function of E)
    e=e.*[1:size(F,2)]; %introduces a face label when edges are counted
    global Fedge; Fedge=zeros(3,size(F,2)); Coordn=[I(:),J(:),e(:)]; %labels in third dimension
    ECNF=ndSparse.build(Coordn(e(:)~=0,:),1,[length(V),length(V),size(F,2)]); %Edge Count on Node, but for each Face (stored in a third dimension)
    Fedge=find3D(findPattern(find3D(ECNF+permute(ECNF,[2,1,3]),'ifo nodes'), E.'),'ifo edges'); %Note: preserves alternating pos/neg orientations
    Fedge([1,2,3],2:2:end)=Fedge([3,1,2],2:2:end); %reshape(E(Fedge(:)),3,[])==F(2,:)+j*F(3,:) %(onschuldige) even permutatie zodat Fedge(1,:)<->F(1,:) en bewijs hiervoor  
    
    %Find boundary edges
    if meshType == 'cylinder'
        %Maar O(N^2) want slechts boundary. 3rde grootste bottleneck na Fedge en eigs
        %TM
        E_vec=round((V(:,imag(E_boundary))-V(:,real(E_boundary)))*(N-1)); %normalized to be logical array. round for floating point error
        E_hori=E_boundary(E_vec(1,:)==1); E_vert=E_boundary(E_vec(2,:)==1);
        E_axis=E_hori(V(2,real(E_hori))==0); E_mantle=E_hori(V(2,real(E_hori))==1); E_leftBound=E_vert(V(1,real(E_vert))==0); E_rightBound=E_vert(V(1,real(E_vert))==1);
        %E_boundary=[E_leftBound,E_rightBound,E_axis,E_mantle];
        %E_boundary_edge=reshape(full(sum(sparse(bsxfun(@eq,E_boundary(:),E.').*[1:length(E)]),2)),[],4);
        E_axis_edge=full(sum(sparse(bsxfun(@eq,E_axis(:),E.').*[1:length(E)]),2));
        E_mantle_edge=full(sum(sparse(bsxfun(@eq,E_mantle(:),E.').*[1:length(E)]),2));
        assocAxisEdge=[];
        for actualAxisEdge=E_axis_edge.'
            assocAxisEdge=[assocAxisEdge;Fedge(:,find(sum(sparse(bsxfun(@eq,Fedge,actualAxisEdge)),1),2))];
        end
        E_axis_edge=assocAxisEdge;
        %E_axis_edge=assocAxisEdge(find(~sum(sparse(bsxfun(@eq,assocAxisEdge.',E_axis_edge)),1)));
        
        %TE
        E_leftBound_node=find(sparse(bsxfun(@eq,V(1,:),0))).';
        E_rightBound_node=find(sparse(bsxfun(@eq,V(1,:),1))).';
        E_axis_node=find(sparse(bsxfun(@eq,V(2,:),0))).';
        E_mantle_node=find(sparse(bsxfun(@eq,V(2,:),1))).';
    end
    
    %RVWn?
    % neumann is altijd goed (kijk heel goed naar interpretatie d/dz e_z
    % met alle drie w_edges en dan op figuur kijken naar pijltjes steeds
    % zelfde z comp, of bewijs wiskundig). Axis is problematisch want daar moet
    % e_r=0 in TM op edges en het beste dat ik kan doen e_z=0 door
    % horizontale edge =0 te doen. Anders kun je zeggen dat het per default al ok is?
    % Dirichlet mantel is te doen door horizontale weg te doen en dan hebt ge e_z=0 zou gewild    
    
    %IDEE: aangezien Er ~ diff(Ez,r) zou je gewoon een neumann RVW kunnen
    %opleggen voor Ez en de boel is opgelost
    %% Additional variables. 
    %Pre-allocate to increase performance
    z=[V(1,F(1,:));V(1,F(2,:));V(1,F(3,:))]; r=[V(2,F(1,:));V(2,F(2,:));V(2,F(3,:))];
    a = zeros(3,size(F,2)); b = zeros(3,size(F,2)); c = zeros(3,size(F,2)); A=zeros(1,size(F,2));
    %matlab can't have function handles in arrays so use cells
    global w_nodes; global w_nodes_total; global w_edges_z; global w_edges_r; global w_edges_z_total; global w_edges_r_total;
    %w_nodes = cell(size(F,2),3); %w_nodes_total = cell(size(F,2),1);
    %w_edges_z = cell(size(F,2),3); w_edges_r = cell(size(F,2),3); 
    %w_edges_z_total = cell(size(F,2),1); w_edges_r_total = cell(size(F,2),1);
    
    i=[1,2,3]; global s; s=@(i) mod(i-1,3)+1; %cyclical indexing

    a(i,:)=z(s(i+1),:).*r(s(i+2),:)-z(s(i+2),:).*r(s(i+1),:);
    b(i,:)=r(s(i+1),:)-r(s(i-1),:);
    c(i,:)=z(s(i-1),:)-z(s(i+1),:);
    A(i,:)=(b(s(i+1),:).*c(s(i+2),:)-b(s(i+2),:).*c(s(i+1),:))/2; %idem for all i
    
    %i in 3rd dimension, n in 4th dimension
    w_nodes=@(z,r,i,n) (permute(a(s(i),n),[3,4,1,2])+permute(b(s(i),n),[3,4,1,2]).*z+permute(c(s(i),n),[3,4,1,2]).*r)./(2*permute(A(i,n),[3,4,1,2]));
    w_edges_z=@(z,r,i,n) (permute(a(s(i+1),n),[3,4,1,2]).*permute(b(s(i+2),n),[3,4,1,2])-permute(a(s(i+2),n),[3,4,1,2]).*permute(b(s(i+1),n),[3,4,1,2])+(permute(c(s(i+1),n),[3,4,1,2]).*permute(b(s(i+2),n),[3,4,1,2])-permute(c(s(i+2),n),[3,4,1,2]).*permute(b(s(i+1),n),[3,4,1,2])).*r)./(4*permute(A(i,n),[3,4,1,2]).^2);
    w_edges_r=@(z,r,i,n) (permute(a(s(i+1),n),[3,4,1,2]).*permute(c(s(i+2),n),[3,4,1,2])-permute(a(s(i+2),n),[3,4,1,2]).*permute(c(s(i+1),n),[3,4,1,2])+(permute(b(s(i+1),n),[3,4,1,2]).*permute(c(s(i+2),n),[3,4,1,2])-permute(b(s(i+2),n),[3,4,1,2]).*permute(c(s(i+1),n),[3,4,1,2])).*z)./(4*permute(A(i,n),[3,4,1,2]).^2);
    
    w_nodes_total=@(z,r,i,n) w_nodes{n,1}(z,r)+w_nodes{n,2}(z,r)+w_nodes{n,3}(z,r);
    w_edges_z_total=@(z,r,i,n) w_edges_z{n,1}(z,r)+w_edges_z{n,2}(z,r)+w_edges_z{n,3}(z,r);
    w_edges_r_total=@(z,r,i,n) w_edges_r{n,1}(z,r)+w_edges_r{n,2}(z,r)+w_edges_r{n,3}(z,r);
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
    %assumpties hier want we weten niet hoe die expansie gaat...
    m_node(1:9,:)=c(i_index,:)/6+c(j_index,:)/6+(c(i_index,:).*c(j_index,:)+b(i_index,:).*b(j_index,:))./(2*abs(A(i_index,:))).*(r(1,:)+r(2,:)+r(3,:))/6;
    M_node=sparse(I,J,m_node,length(V),length(V));
        
    %BCs    
    sz1=size(G_edge,1); sz2=size(E_mantle_edge,1);
    G_edge(E_mantle_edge,:)=zeros(sz2,sz1);
    G_edge(:,E_mantle_edge)=zeros(sz1,sz2);
    M_edge(E_mantle_edge,:)=zeros(sz2,sz1);
    M_edge(:,E_mantle_edge)=zeros(sz1,sz2);
    for n=1:size(E_mantle_edge,1)
        M_edge(E_mantle_edge(n),E_mantle_edge(n))=1;
        %G_edge(E_mantle_edge(n),E_mantle_edge(n))=1;
    end
    
    sz1=size(G_edge,1); sz2=size(E_axis_edge,1);
    G_edge(E_axis_edge,:)=zeros(sz2,sz1);
    G_edge(:,E_axis_edge)=zeros(sz1,sz2);
    M_edge(E_axis_edge,:)=zeros(sz2,sz1);
    M_edge(:,E_axis_edge)=zeros(sz1,sz2);
    for n=1:size(E_axis_edge,1)
        M_edge(E_axis_edge(n),E_axis_edge(n))=1;
        %G_edge(E_axis_edge(n),E_axis_edge(n))=1;
    end

    sz1=size(G_node,1); sz2=size(E_leftBound_node,1);
    G_node(E_leftBound_node,:)=zeros(sz2,sz1);
    G_node(:,E_leftBound_node)=zeros(sz1,sz2);
    M_node(E_leftBound_node,:)=zeros(sz2,sz1);
    M_node(:,E_leftBound_node)=zeros(sz1,sz2);
    for n=1:size(E_leftBound_node,1)
        M_node(E_leftBound_node(n),E_leftBound_node(n))=1;
    end
    
    sz1=size(G_node,1); sz2=size(E_rightBound_node,1);
    G_node(E_rightBound_node,:)=zeros(sz2,sz1);
    G_node(:,E_rightBound_node)=zeros(sz1,sz2);
    M_node(E_rightBound_node,:)=zeros(sz2,sz1);
    M_node(:,E_rightBound_node)=zeros(sz1,sz2);
    for n=1:size(E_rightBound_node,1)
        M_node(E_rightBound_node(n),E_rightBound_node(n))=1;
    end
    
    sz1=size(G_node,1); sz2=size(E_axis_node,1);
    G_node(E_axis_node,:)=zeros(sz2,sz1);
    G_node(:,E_axis_node)=zeros(sz1,sz2);
    M_node(E_axis_node,:)=zeros(sz2,sz1);
    M_node(:,E_axis_node)=zeros(sz1,sz2);
    for n=1:size(E_axis_node,1)
        M_node(E_axis_node(n),E_axis_node(n))=1;
    end
    
    sz1=size(G_node,1); sz2=size(E_mantle_node,1);
    G_node(E_mantle_node,:)=zeros(sz2,sz1);
    G_node(:,E_mantle_node)=zeros(sz1,sz2);
    M_node(E_mantle_node,:)=zeros(sz2,sz1);
    M_node(:,E_mantle_node)=zeros(sz1,sz2);
    for n=1:size(E_mantle_node,1)
        M_node(E_mantle_node(n),E_mantle_node(n))=1;
    end

    %length(E)-size(F,2)+1 %eig gwn #nodes volgens euler
    %svd(full(M_edge)); %sum(svd(full(M_edge))<eps('single'))
    %sum(svd(full(M_edge))<eps('single') | abs(svd(full(M_edge))-1)<eps('single'))
    %[Eedge,D]=eig(full(M_edge),500*full(G_edge)); %blijkbaar zou dit een juister resultaat geven dan eigs
    %myeigs=eig(full(M_edge),full(G_edge));
    %sum(myeigs<1)
    %plot([1:1:size(sort(myeigs))],sort(myeigs))
    %global TM_modes; [TM_modes,D_TM]=eigs(M_edge,G_edge,5,24.11); 
    %global TE_modes; [TE_modes,D_TE]=eigs(M_node,G_node,5,5);
    global TM_modes; [TM_modes,D_TM]=eig(full(M_edge),full(G_edge));
    global TE_modes; [TE_modes,D_TE]=eig(full(M_node),full(G_node));
    
    %sort the eigenvalues and associated eigenmodes from low to high 
    [dnnz_TM,sortingIndices]=sort(D_TM(D_TM>eps('single'))); TM_modes=TM_modes(:,sortingIndices);
    [dnnz_TE,sortingIndices]=sort(D_TE(D_TE>eps('single'))); TE_modes=TE_modes(:,sortingIndices);
    dnnz_TM
    dnnz_TE
    
    disp(['TM resonant frequency :' num2str(dnnz_TM(mode))])
    disp(['TE resonant frequency :' num2str(dnnz_TE(mode))])
    
    %%{
    if meshType == 'cylinder'
        if (Nz+Nr)/2 < 61
            res=1/(N-1)/6; %roughly 6 per triangle length. Increasing res betters low N visualization but not eigenvalues
            %res=0.02; %best to pick a round number since the edges and nodes are located on "round" places
            viewElectricField(mode, res, z_max, r_max)
            %viewMeshandBasisFcts(1, 1, res,z_max,r_max)
            %viewMeshandBasisFcts(2, 1, res,z_max,r_max)
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

function out=trianglesImIn(Z,R) 
    global F; global w_nodes; 
    %when visualizing, this is a clear bottleneck as the current algorithm
    %uses an expensive point-in-poly algorithm. One could make this
    %O(n*log(n)) with e.g. binary search but that seems quite complicated.
    %The current algorithm visualizes N=50 in 40s and we know the real limits
    %are around N=150 with eigs and find3D routine. But what really is the
    %point of visualization at this point? Errors? No, we only do that wrt
    %the eigenvalues
    
    %vectorized algorithm but sadly takes too much memory
    %data=w_nodes(Z,R,1:3,1:size(F,2)); out=ndSparse(data(:,:,1,:)<=1+eps('single') & data(:,:,1,:)>=0-eps('single') & data(:,:,2,:)<=1+eps('single') & data(:,:,2,:)>=0-eps('single') & data(:,:,3,:)<=1+eps('single') & data(:,:,3,:)>=0-eps('single'));
    
    %so that's why I have to do a for loop...
    %{
    out=ndSparse.spalloc([size(Z), size(F,2)],size(Z,1)^2);
    for n=1:size(F,2)
        data=w_nodes(Z,R,1:3,n);
        out(:,:,n)=sparse(data(:,:,1)<=1+eps('single') & data(:,:,1)>=0-eps('single') & data(:,:,2)<=1+eps('single') & data(:,:,2)>=0-eps('single') & data(:,:,3)<=1+eps('single') & data(:,:,3)>=0-eps('single'));
    end
    %}

    %mixed form with chunks
    availableMemory=1.5*10^9; availableElements=availableMemory/8; %1 element = 1 byte and my pc doesn't have a lot of RAM
    requestedElements=size(Z,1)^2*3*size(F,2); amountOfChunks=requestedElements/availableElements; chunkSize=floor(size(F,2)/amountOfChunks);
    if requestedElements < availableElements
        data=w_nodes(Z,R,1:3,1:size(F,2)); out=ndSparse(data(:,:,1,:)<=1+eps('single') & data(:,:,1,:)>=0-eps('single') & data(:,:,2,:)<=1+eps('single') & data(:,:,2,:)>=0-eps('single') & data(:,:,3,:)<=1+eps('single') & data(:,:,3,:)>=0-eps('single'));
    else
        out=ndSparse.spalloc([size(Z), size(F,2)],size(Z,1)^2);
        disp(['Initiating visualization processing with ' num2str(floor(amountOfChunks)+1) ' ' num2str(availableMemory*10^(-9)) ' GB data chunks'])
        for n=1:floor(amountOfChunks)
            disp(['Starting chunk ' num2str(n) ' of ' num2str(floor(amountOfChunks)+1)])
            chunk=((n-1)*chunkSize+1):(n*chunkSize);
            data=w_nodes(Z,R,1:3,chunk); out(:,:,chunk)=ndSparse(data(:,:,1,:)<=1+eps('single') & data(:,:,1,:)>=0-eps('single') & data(:,:,2,:)<=1+eps('single') & data(:,:,2,:)>=0-eps('single') & data(:,:,3,:)<=1+eps('single') & data(:,:,3,:)>=0-eps('single'));
        end
        disp(['Starting chunk ' num2str(floor(amountOfChunks)+1) ' of ' num2str(floor(amountOfChunks)+1)])
        finalChunk=(floor(amountOfChunks)*chunkSize+1):size(F,2);
        data=w_nodes(Z,R,1:3,finalChunk); out(:,:,finalChunk)=ndSparse(data(:,:,1,:)<=1+eps('single') & data(:,:,1,:)>=0-eps('single') & data(:,:,2,:)<=1+eps('single') & data(:,:,2,:)>=0-eps('single') & data(:,:,3,:)<=1+eps('single') & data(:,:,3,:)>=0-eps('single'));
    end
    %sadly this doesn't seem to give significant improvement. I don't know
    %why
end

function Ephi=ElectricField_TE(mode,Z,R,triList)
    global w_nodes; global F; global TE_modes
    Ephi=zeros(size(Z,1),size(R,2));
    
    for n=1:size(F,2) %we only evaluate the function on local support so operation is O(n*1/n) cost effective        
        nnzZRindex=find(triList(:,:,n)); nnzZCoord=Z(nnzZRindex); nnzRCoord=R(nnzZRindex);
        Ephi(nnzZRindex)=sum(reshape(TE_modes(F(:,n),mode),1,1,3,1).*w_nodes(nnzZCoord,nnzRCoord,1:3,n),3); %don't average as it's perfectly continuous
    end
end

function [Ez,Er]=ElectricField_TM(mode,Z,R,triList)
    global w_edges_z; global w_edges_r; global Fedge; global F; global TM_modes
    Ez=zeros(size(Z)); Er=zeros(size(Z));
    
    triAmounts=sum(triList,4); %for averaging effect on edges    
    for n=1:size(F,2) %we only evaluate the function on local support so operation is O(n*1/n) cost effective        
        nnzZRindex=find(triList(:,:,n)); nnzZCoord=Z(nnzZRindex); nnzRCoord=R(nnzZRindex);
        Ez(nnzZRindex)=Ez(nnzZRindex)+sum(reshape(TM_modes(Fedge(:,n),mode),1,1,3,1).*w_edges_z(nnzZCoord,nnzRCoord,1:3,n),3)./triAmounts(nnzZRindex);
        Er(nnzZRindex)=Er(nnzZRindex)+sum(reshape(TM_modes(Fedge(:,n),mode),1,1,3,1).*w_edges_r(nnzZCoord,nnzRCoord,1:3,n),3)./triAmounts(nnzZRindex);
    end
end

function out = find3D(X,ifo) 
    [~,y,z]=size(X);
    flat=mod(find(permute(ndSparse(X),[2,1,3]))-1,y)+1; %let op: doet al een triu (idk hoe?)
    if ifo == 'ifo nodes'
        out=reshape(flat(1:2:end)+1i*flat(2:2:end),[3,1,z]); %3d shape
    elseif ifo == 'ifo edges'
        out=reshape(flat,[3,z]); %2d shape
    end
end

function out = find2D(X) %for elegance in one line and efficient double-element storage
    [x,y]=find(X); out = (x+1i*y);
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

function viewElectricField(mode,resolution,z_max,r_max)
    global V; global F; global N;
    [Z,R] = meshgrid(0:resolution:z_max,0:resolution:r_max); triList=trianglesImIn(Z,R);
    [Ez,Er]=ElectricField_TM(mode,Z,R,triList); %Eabs=sqrt(Ez.^2+Er.^2);
    Ephi=ElectricField_TE(mode,Z,R,triList);    
    
    figure('Name',['Ez in TM mode ' num2str(mode) ' (N = ' num2str(N) ')'])
    hold on
    if N <= 30
        V=V.'; F=F.'; patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F,'FaceColor','none','EdgeColor','black','LineWidth',1.5); V=V.'; F=F.';     
    end
    surface(Z,R,Ez); view(30,90-30); 
    %quiver(Z,R,Ez,Er,'black'); 
    alpha 0.75
    xlabel('$z$ [$m$]','interpreter','latex')
    ylabel('$r$ [$m$]','interpreter','latex')
    hold off

    figure('Name',['Er in TM mode ' num2str(mode) ' (N = ' num2str(N) ')'])
    hold on
    if N <= 30
        V=V.'; F=F.'; patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F,'FaceColor','none','EdgeColor','black','LineWidth',1.5); V=V.'; F=F.';     
    end
    surface(Z,R,Er); view(30,90-30); 
    %quiver(Z,R,Ez,Er,'black'); 
    alpha 0.75
    xlabel('$z$ [$m$]','interpreter','latex')
    ylabel('$r$ [$m$]','interpreter','latex')
    hold off
    
    figure('Name',['Ephi in TE mode ' num2str(mode) ' (N = ' num2str(N) ')'])
    hold on
    if N <= 30
        V=V.'; F=F.'; patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F,'FaceColor','none','EdgeColor','black','LineWidth',1.5); V=V.'; F=F.';     
    end
    surface(Z,R,Ephi); view(30,90-30); 
    %quiver(Z,R,Ephi,Ephi,'black'); 
    alpha 0.75
    xlabel('$z$ [$m$]','interpreter','latex')
    ylabel('$r$ [$m$]','interpreter','latex')
    hold off
end

function viewMeshandBasisFcts(triangle,point,resolution,z_max,r_max)
    global V; global F; global w_nodes; global w_nodes_total; global w_edges_z; global w_edges_r; global w_edges_z_total; global w_edges_r_total; global N
    V=V.'; %to do: write viewMesh so this isn't necessary(?)
    F=F.';
    [Z,R] = meshgrid(0:resolution:z_max,0:resolution:r_max);

    figure('Name',['Mesh and basis function(s) of triangle ' num2str(triangle) ' and node ' num2str(point)  ' (N = ' num2str(N) ')'])
    hold on
    if point ~= 'total'
        %show w_node as filled contourplot
        [c,h] = contourf(Z,R,w_nodes{triangle,point}(Z,R),1/(2*resolution),'--','ShowText','on');        
        h.LevelList=round(h.LevelList,2);  %rounds levels to 2nd decimal place
        clabel(c,h);
        %show w_edge as vectorplot
        quiver(Z,R,w_edges_z{triangle,point}(Z,R),w_edges_r{triangle,point}(Z,R),'black');
    else
        point=1;
        [c,h] = contourf(Z,R,round(w_nodes_total{triangle,point}(Z,R),3),1/(2*resolution),'--','ShowText','on');                
        h.LevelList=round(h.LevelList,2);  %rounds levels to 2nd decimal place
        clabel(c,h);
        quiver(Z,R,w_edges_z_total{triangle,point}(Z,R),w_edges_r_total{triangle,point}(Z,R),'black'); 
    end
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