function areas16(a,dx,dy,n)
tol = 0.000000000000000001;
global mx16 my16 area locx locy
area = [];
cc = n/2+0.5;
for v=1:size(mx16,2)
    i = mx16(1,v);
    j = my16(1,v);
    x = (i-cc)*dx;
    y = (cc-j)*dy;
    x1=x-dx/2;
    x2=x+dx/2;
    y1=y-dy/2;
    y2=y+dy/2;
    s = dx/2000;
    lx = [x1+s/2:s:x2-s/2];
    ly = [y1+s/2:s:y2-s/2];
    tt = dx/s;
    opp = 0;
        for ii=1:tt
            for jj=1:tt
                if (lx(ii)-locx(1,v))^2+(ly(jj)-locy(1,v))^2 < a^2
                    opp=opp+s^2;
                end
            end
        end
        opp = dx*dy-opp;
    if opp <= 0.005*dx*dy
        opp = 0;
    end
    area = [area opp];
end