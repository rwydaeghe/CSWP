function positions(a,dx,dy,n)
tol = 0.000000000000001;
cc = n/2+0.5;
global mx my
mx=[];
my=[];
for i = 1:n
    for j = 1:n
        x = (i-cc)*dx;
        y = (cc-j)*dy;
        x1=x-dx/2;
        x2=x+dx/2;
        y1=y-dy/2;
        y2=y+dy/2;
        if x1^2+y1^2 <= a^2 || x2^2+y1^2 <= a^2 || x1^2+y2^2 <= a^2 || x2^2+y2^2 <= a^2
            if x1^2+y1^2 >= a^2 & x2^2+y1^2 >= a^2 & x1^2+y2^2 >= a^2 & x2^2+y2^2 >= a^2
            else
                mx = [mx i];
                my = [my j];
            end
        end
    end
end



                    
            