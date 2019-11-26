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