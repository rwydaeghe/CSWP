function viewMeshandBasisFcts(triangle,point,resolution,z_max,r_max)
    global V; global F; global w_nodes; global w_edges_z; global w_edges_r;
    [Z,R] = meshgrid(0:resolution:z_max,0:resolution:r_max);
    
    hold on
    %show w_node as filled contourplot
    [c,h] = contourf(Z,R,w_nodes{triangle,point}(Z,R),1/(2*resolution),'--','ShowText','on');
    h.LevelList=round(h.LevelList,2);  %rounds levels to 2nd decimal place
    clabel(c,h);
    %show w_edge as vectorplot
    %quiver(Z,R,w_edges_z{triangle,point}(Z,R),w_edges_r{triangle,point}(Z,R),'black'); 
    %Z,R
    %Ez(Z,R)
    quiver(Z,R,Ez(Z,R),Er(Z,R),'black'); 
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

