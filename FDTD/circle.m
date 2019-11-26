function circle(dx,dy,a) %n even
global ox oy p 
n = size(p,1);
% n=100;
m = n*dx/2;
% q = ones(100,100)
if mod(n,2) == 0
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
                p(j,i) = 0;
                ox(j,i) = 0;
                ox(j,i+1) = 0;
                oy(j,i) = 0;
                oy(j+1,i) = 0;
            end
        end
    end
end
%gridplot(q)
end

