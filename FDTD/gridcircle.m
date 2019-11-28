dx=0.5;
dy=dx;
a=2;
ri=2;
ro=2+4*dx;
dh=dx/ro;
r=[ri:dx:ro];
h=[0:dh:2*pi];
f=[0:0.001:2*pi];
z = [-2.5*a:dx:2.5*a];
n=5*a/dx;
q=ones(n);
c = n/2+0.5;
for i = 1:n
        for j = 1:n
            x = (i-c)*dx;
            y = (c-j)*dy;
            if x^2+y^2 < ri^2
                q(j,i)=0;
            elseif x==0
                l=min([y-ri y-(-ri)]);
                if l <= 1.25*dx
                    q(j,i)=0;
                end
            else
                xi=sqrt(ri^2/(1+(y/x)^2));
                yi=y/x*xi;
                l=min([sqrt((y-yi)^2+(x-xi)^2) sqrt((y+yi)^2+(x+xi)^2)]);
                if l <= 1.25*dx
                    q(j,i)=0;
                end
            end
                
        end
end



figure (1);
hold on
axis equal
for i = 1:length(r)
    plot(r(i)*cos(f),r(i)*sin(f),'k');
end

for i = 1:length(h)
    plot([ri ro]*cos(h(i)),[ri ro]*sin(h(i)),'k');
end
for i = 1:length(z)
         plot(z, z(i)*ones(length(z),1),'k');
         plot(z, z(i)*ones(length(z),1),'k');
          plot(z(i)*ones(length(z),1), z,'k');
end

hold off

figure (2)
gridplot(q)



    



