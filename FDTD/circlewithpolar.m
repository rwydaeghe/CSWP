function circlewithpolar(dx,dy,a,n)
global q qx qy
ri = a;
c = n/2+0.5;
for i = 1:n
    for j = 1:n
        x = (i-c)*dx;
        y = (c-j)*dy;
        if x^2+y^2 < (ri+1.25*dx)^2
                q(j,i)=0;
                qx(j,i) = 0;
                qx(j,i+1) = 0;
                qy(j,i) = 0;
                qy(j+1,i) = 0;
        end
    end
end