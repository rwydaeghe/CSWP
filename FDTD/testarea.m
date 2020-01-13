function areas(a,dx,dy,n)
global mx my area
area = [];
for v=1:size(mx,2)
    i = mx(1,v)
    j = my(1,v)
    x = (i-c)*dx;
    y = (c-j)*dy;
    x1=x-dx/2;
    x2=x+dx/2;
    y1=y-dy/2;
    y2=y+dy/2;
    s = dx/2000;
    lx = [x1:s:x2];
    ly = [y1:s:y2];
    tt = dx/s+1;
    opp = 0;
    if x1^2+y1^2 <= a^2 & x2^2+y1^2 <= a^2 & x1^2+y2^2 <= a^2 & x2^2+y2^2 <= a^2
    else
        for ii=1:tt
            for jj=1:tt
                if lx(ii)^2+ly(jj)^2 <=a^2
                    opp=opp+s^2;
                end
            end
        end
        opp = dx*dy-opp;
    end
    area = [area opp];
end