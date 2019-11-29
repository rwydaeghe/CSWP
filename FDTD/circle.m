function circle(dx,dy,a,n) %n even
global q qx qy
q = ones(n,n);
qx = ones(n,n+1);
qy = ones(n+1,n);
c = n/2+0.5;
for i = 1:n
    for j = 1:n
        x = (i-c)*dx;
        y = (c-j)*dy;
%             if x > 0
%                 x = x-dx/2;
%             else
%                 x = x+dx/2;
%             end
%             if y > 0
%                 y = y-dy/2;
%             else
%                 y = y+dy/2;
%             end
        if x^2+y^2 < a^2
%                 q(j,i)=0;
            q(j,i) = 0;
            qx(j,i) = 0;
            qx(j,i+1) = 0;
            qy(j,i) = 0;
            qy(j+1,i) = 0;
        end
    end
end
%gridplot(q)
end

