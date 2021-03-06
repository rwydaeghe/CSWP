%%% List of kinda usefull code %%%

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
%F in fct van E:
    %alternatief, niet gevectorizeerd:
        %Fedge(1,u)=(strfind(reshape(E,1,[]),[x(1),y(1)])+1)/2;
        %Fedge(2,u)=(strfind(reshape(E,1,[]),[x(2),y(2)])+1)/2;
        %Fedge(3,u)=(strfind(reshape(E,1,[]),[x(3),y(3)])+1)/2;
%thanks to https://nl.mathworks.com/matlabcentral/answers/317880-how-to-vectorize-strfind
%S = [1+j*1,1+j*2,2+3*j,1,5+j*4,0.5+3*j,j*1]
%targets = [1+j*2; j*1];
%output = [targets(:,1) == S]
%[~,y]=find(output); y
%i=[1,2,3];
%n=10;
%w=zeros(n,3)
%w(:,1)='hi'


%F edges attempts and successes:

    
%attempt ndSparse vectorially O(n)
%full(speye(length(F))).*ones(length(F),length(F),9)
%working mask:
%Coord=[zeros(length(F),1)+2,ones(length(F),2).*[1:length(F)].';
%       zeros(length(F),1)+3,ones(length(F),2).*[1:length(F)].';
%       zeros(length(F),1)+6,ones(length(F),2).*[1:length(F)].'];
%masked_e=ndSparse.build(Coord,1,[9,length(F),length(F)]);

%I=I.*ones(9,length(F),length(F))
%J=J.*ones(9,length(F),length(F))
%e=reshape([zeros(length(F)),speye(length(F)),speye(length(F)),zeros(length(F)),zeros(length(F)),speye(length(F)),zeros(length(F)),zeros(length(F)),zeros(length(F))],[9,length(F),length(F)]);
%full(e)
%e=ndSparse.build([:;:;2,3,6],speye(length(F)),[length(F),length(F),9])
%e=zeros(9,length(F),length(F)); e([2,3,6],:,:)=1; 
%ME2=ndSparse.build(I,J,e.*mask(u,:),length(V),length(V));
%OLD O(n^2) algo:
%e=zeros(9,length(F)); e([2,3,6],:)=1;
%for u=1:length(F)
%    ME2=sparse(I,J,e.*mask(u,:),length(V),length(V));
%    Fedge(:,u) = imag(find2D([find2D(triu(ME2+ME2.'))==E.']));
%end
%for u = 1:length(F)
%    Fedge = imag(find2D([find2D(L)==E.']))
%    Fedge(:,u) = imag(find2D([find3D(ME2+permute(ME2,[2,1,3]))==E.']));
%    Fedge(:,u) = imag(find2D(bsxfun(@(x,y) x==y, find3D(ME2+permute(ME2,[2,1,3])), E.')))
%end

    %Analyse visualisatie voor windows:
    %{
    hold on
    plot([1:length(Fedge)],mean(Fedge))
    hold on
    plot([1:length(Fedge)],max(Fedge)-min(Fedge))
    m=mean(Fedge); std=max(Fedge)-min(Fedge); 
    output=[m(end)/length(F),max(std)/sqrt(length(F))];
    %Conclusions:
    %the evolution of the means goes as 1:2.5*length(F) with F
    %the max of (max-min) is <2*sqrt(length(F))
    %}



%findPattern
%Bcell=zeros(length(F),3);
    %for u=1:length(F)
    %    u
    %    Bcell(u,:)
    %    Bcell(u,:)=mat2cell(B,1,[windowstart(u)-1,windowsize(u),length(B)-windowend(u)+1])
    %    %inter(u,:)=Bcell{u,2};
    %    %Bcell{u,3}=mat2cell(B,1,[windowstart(u),windowsize(u),length(B)-windowend(u)])
    %end
    %r=Bcell(1,:)
    %m=r{1,1}
    %m{1,2}
    
%B=B.*ones(1,length(B),lenF);
%B2=zeros(1,length(B),lenF);
%out_prov=zeros(3,length(B),lenF); %provisoir                

%in for loop

%B2(1,window1start(u):window1end(u),u)=B(1,window1start(u):window1end(u));
%B2(1,window2start(u):window2end(u),u)=B(1,window2start(u):window2end(u));
%B(1,1:window1start(u)-1,u)=0; B(1,window1end(u)+1:window2start(u)-1,u)=0; B(1,window2end(u)+1:end,u)=0;
%out_prov(:,:,u)=sparse(bsxfun(@eq,A(:,:,u),B(:,:,u)));


%used as an alternative to visualize the electric field without weights
%{
if size(tri,2) == 1
    for edge=1:3
        out(iz,ir)=out(iz,ir)+w_edges_z{tri,edge}(z,r);
    end
elseif size(tri,2) == 2
    commonGlobalEdge=intersect(Fedge(:,tri(1)),Fedge(:,tri(2))); notThis1=real(E(commonGlobalEdge)); notThis2=imag(E(commonGlobalEdge));
    firstTriangleOppositeNode=find((F(:,tri(1))~=notThis1 & F(:,tri(1))~=notThis2));
    secondTriangleOppositeNode=find((F(:,tri(2))~=notThis1 & F(:,tri(2))~=notThis2));
    out(iz,ir)=w_edges_z{tri(1),firstTriangleOppositeNode}(z,r)+w_edges_z{tri(2),secondTriangleOppositeNode}(z,r);
end
%}

%poging tot vectorizeren van buitenste for loop als je Ez,Er,Ephi
%berekent... heeft geen nut gehad want kost met for loop is toch al
%O(n*1/n)

%omslachtiger perpAxisEdge berekenen
        %{
        assocAxisEdge=[]; %all three edges associated with the triangle that contains the E_axis_edge 
        for actualAxisEdge=E_axis_edge.'
            assocAxisEdge=[assocAxisEdge;Fedge(:,find(sum(sparse(bsxfun(@eq,Fedge,actualAxisEdge)),1),2))]
        end
        perpAxisEdge=[];
        for assocEdge=assocAxisEdge.'
            assocEdgeVec=round((V(:,imag(E(assocEdge)))-V(:,real(E(assocEdge))))*(N-1));
            if assocEdgeVec(2)==1 
                perpAxisEdge=[perpAxisEdge;assocEdge];
            end
        end
        %}