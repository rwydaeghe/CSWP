function viewMesh(V,F,varargin)
V = [V, zeros(size(V,1), 1)];
patch('Vertices',V,'Faces',F,'FaceColor','red','EdgeColor','black',varargin{:})
xlabel('x')
ylabel('y')
return