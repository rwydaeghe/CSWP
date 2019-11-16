function circle(dx,dy,a) %n even
global ox oy p 
n = size(p,1);
m = n*dx/2;
if mod(n,2) == 0
    c = n/2+0.5;
    for i = 1:n
        for j = 1:n
            x = (j-c)*dx;
            y = (c-i)*dy;
            if x > 0
                x = x-dx/2;
            else
                x = x+dx/2;
            end
            if y > 0
                y = y-dy/2;
            else
                y = y+dy/2;
            end
            if x^2+y^2 < a^2
                p(i,j) = 0;
                ox(i,j) = 0;
                ox(i,j+1) = 0;
                oy(i,j) = 0;
                oy(i,j+1) = 0;
            end
        end
    end
end
end

