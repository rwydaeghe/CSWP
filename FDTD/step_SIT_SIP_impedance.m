function step_SIT_SIP_impedance(nx,ny,c,dx,dy,dt,a,nc)
global ox oy p
global qx qy q
global oxl oyl plo plp
global oxr oyr pro prp
global oxu oyu puo pup
global oxd oyd pdo pdp
Z=c; %surface impedance for non-reflecting boundary

km = 10000;
n=100;

for i=1:nx
    for j=1:n
        k = km*(((j-1/2)*dy)/n/dy)^4;
        if j == 1
            puo(j,i) = ((1-k/2*dt)*puo(j,i) - c^2*dt/dy*(oyu(j,i)-oy(ny+1,i)))/(1+k/2*dt);
        else
            puo(j,i) = ((1-k/2*dt)*puo(j,i) - c^2*dt/dy*(oyu(j,i)-oyu(j-1,i)))/(1+k/2*dt);
        end
        pup(j,i) = pup(j,i) - c^2*dt/dx*(oxu(j,i+1)-oxu(j,i));
    end
end

for i=1:nx
    for j=1:ny
        p(j,i) = p(j,i) - c^2*dt*(1/dx*(ox(j,i+1)-ox(j,i))+1/dy*(oy(j+1,i)-oy(j,i)));
    end
end

for i=2:nx
    for j=1:n
        oxu(j,i) = oxu(j,i) - dt*(1/dx*(pup(j,i)+puo(j,i)-pup(j,i-1)-puo(j,i-1)));
    end
end

for j=1:n
    oxu(j,1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*oxu(j,1)-2*dt/dx*(pup(j,1)+puo(j,1)));
    oxu(j,nx+1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*oxu(j,nx+1)+2*dt/dx*(pup(j,nx)+puo(j,nx)));
end

for i=1:nx
    for j=1:n
        k = km*(((j)*dy)/n/dy)^4;
        if j == n
            oyu(n,i) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oyu(n,i)+2*dt/dx*(pup(n,i)+puo(n,i))); %look
        else
            oyu(j,i) = ((1-k/2*dt)*oyu(j,i) - dt*(1/dy*(pup(j+1,i)+puo(j+1,i)-pup(j,i)-puo(j,i))))/(1+k/2*dt);
        end
    end
end


for i=2:nx
    for j=1:ny
        ox(j,i) = ox(j,i) - dt*(1/dx*(p(j,i)-p(j,i-1))); 
    end
end

for j=1:ny
    ox(j,1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(j,1)-2*dt/dx*p(j,1));
    ox(j,nx+1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(j,nx+1)+2*dt/dx*p(j,nx));
end

for i=1:nx
    for j = 2:ny
        oy(j,i) = oy(j,i) - dt*(1/dy*(p(j,i)-p(j-1,i)));
    end
end

for i=1:nx
    oy(1,i) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(1,i)-2*dt/dx*p(1,i));
    oy(ny+1,i) = oy(ny+1,i) - dt/dy*(pup(1,i)+puo(1,i)-p(ny,1));
end
% oy = oy.*qy;
% ox = ox.*qx;
% p = p.*q;


end
        








%update ox oy p

 