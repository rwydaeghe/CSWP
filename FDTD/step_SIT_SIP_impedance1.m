function step_SIT_SIP_impedance1(nx,ny,c,dx,dy,dt)
global ox oy p
Z=c; %surface impedance for non-reflecting boundary

for i=1:nx
    for j=1:ny
        p(i,j) = p(i,j) - c^2*dt*(1/dx*(ox(i+1,j)-ox(i,j))+1/dy*(oy(i,j+1)-oy(i,j)));
    end
end

for i=2:nx
    for j=2:ny
        ox(i,j) = ox(i,j) - dt*(1/dx*(p(i,j)-p(i-1,j)));
        oy(i,j) = oy(i,j) - dt*(1/dy*(p(i,j)-p(i,j-1)));
    end
end
for j=2:ny
    ox(1,j) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(1,j)-2*dt/dx*p(1,j));
    ox(nx+1,j) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(nx+1,j)+2*dt/dx*p(nx,j));
end
for i=2:nx
    oy(i,1) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(i,1)-2*dt/dx*p(i,1));
    oy(i,ny+1) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(i,ny+1)+2*dt/dx*p(i,ny));
end


end
        








%update ox oy p

 