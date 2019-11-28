function step_SIT_SIP_impedance(nx,ny,c,dx,dy,dt,a,nc)
global ox oy p
Z=c; %surface impedance for non-reflecting boundary

for i=1:nx
    for j=1:ny
        p(j,i) = p(j,i) - c^2*dt*(1/dx*(ox(j,i+1)-ox(j,i))+1/dy*(oy(j+1,i)-oy(j,i)));
    end
end

for i=2:nx
    for j=2:ny
        ox(j,i) = ox(j,i) - dt*(1/dx*(p(j,i)-p(j,i-1)));
        oy(j,i) = oy(j,i) - dt*(1/dy*(p(j,i)-p(j-1,i)));
    end
end
for j=2:ny
    ox(j,1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(j,1)-2*dt/dx*p(j,1));
    ox(j,nx+1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(j,nx+1)+2*dt/dx*p(j,nx));
end
for i=2:nx
    oy(1,i) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(1,i)-2*dt/dx*p(1,i));
    oy(ny+1,i) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(ny+1,i)+2*dt/dx*p(ny,i));
end
circle(dx,dy,a)


end
        








%update ox oy p

 