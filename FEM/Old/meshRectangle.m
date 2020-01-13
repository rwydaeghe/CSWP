function [vertices, faces]=meshRectangle(lower_left, upper_right, Nx, Ny)

% lower_left en upper_right are 2D vectors with coordinates corresponding to the
% corners of the rectang;e.
% Nx: number of nodes along the horizontal direction
% Ny: number of nodes along the vertical direction
% vectices: an Nx*Ny x 2 matrix containing the coordinates of the nodes
% faces: a 2*(Nx-1)*(Ny-1) x 3 matrix containing the nodes per face

delta_x = (upper_right(1)-lower_left(1))/(Nx-1);
delta_y = (upper_right(2)-lower_left(2))/(Ny-1);

%Nx = (upper_right(1)-lower_left(1))/delta_x + 1;
%Ny = (upper_right(2)-lower_left(2))/delta_y + 1;

Xvector = [lower_left(1):delta_x:upper_right(1)];
Yvector = [lower_left(2):delta_y:upper_right(2)];

vertices = zeros(Nx*Ny, 3);
for p=1:Nx
    for q=1:Ny
        vertices(q+(p-1)*Ny,1)=Xvector(p);
        vertices(q+(p-1)*Ny,2)=Yvector(q);
        vertices(q+(p-1)*Ny,3)=0;
    end
end

faces = zeros(2*(Nx-1)*(Ny-1), 3);
index = 1;
for p=1:Nx-1
    for q=1:Ny-1
        faces(2*index-1,1) = q+(p-1)*Ny;
        faces(2*index-1,2) = q+p*Ny;
        faces(2*index-1,3) = q+1+(p-1)*Ny;
        index = index + 1;
    end
end
index = 1;
for p=2:Nx
    for q=2:Ny
        faces(2*index,1) = q+(p-1)*Ny;
        faces(2*index,2) = q+(p-2)*Ny;
        faces(2*index,3) = q-1+(p-1)*Ny;
        index = index + 1;
    end
end


return