function viewMeshandBasisFcts(V,F,w_nodes,w_edges,triangle,point,resolution)
    [Z,R] = meshgrid(0:resolution:1,0:resolution:1);
    
    %unpacking w_edges onto the mesh needs special attention
    wz = zeros(length(Z),length(R));
    wr = zeros(length(Z),length(R));
    for i = 1:length(Z)
        for j = 1:length(R)
            wzr = w_edges{triangle,point}(Z(i),R(j)); %gedoe moet want anders verwart matlab indexing met functie-argumenten
            wz(i,j) = wzr(1);
            wr(i,j) = wzr(2);            
        end
    end
    
    hold on
    %show w_node as filled contourplot
    [c,h] = contourf(Z,R,w_nodes{triangle,point}(Z,R),1/(2*resolution),'--','ShowText','on');
    h.LevelList=round(h.LevelList,2);  %rounds levels to 2nd decimal place
    clabel(c,h);
    %show w_edge as vectorplot
    quiver(Z,R,wz,wr,'black');
    %show mesh
    patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F,'FaceColor','none','EdgeColor','black','LineWidth',1.5);
    %highlight selected triangle in blue
    patch('Vertices',[V, zeros(size(V,1), 1)],'Faces',F(triangle,:),'FaceColor','none','EdgeColor','blue', 'LineWidth', 1.5);
    %highlight selected node in red
    scatter(V(F(triangle,point),1), V(F(triangle,point),2),'filled','red');        
    %highlight selected edge in red
    point = point - 1; %indices in matlab sucken weer
    p2=mod(point+1,3)+1;
    p3=mod(point+2,3)+1;
    plot([V(F(triangle,p2),1), V(F(triangle,p3),1)],[V(F(triangle,p2),2), V(F(triangle,p3),2)],'red','LineWidth',1.5);
    xlabel('$z$','interpreter','latex')
    ylabel('$r$','interpreter','latex')
    hold off
end

