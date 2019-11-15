function viewMesh(V,F,varargin)
V = [V, zeros(size(V,1), 1)];
V = V(:,1:2); %moet er blijkbaar bij? 
patch('Vertices',V,'Faces',F,'FaceColor','red','EdgeColor','black',varargin{:})
xlabel('x')
ylabel('y')
return