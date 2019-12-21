function interpolouteredgecircle(Nr,Nh,ro,dr,dh,n,dx,dy)
global or
global ox oy
c = n/2+0.5;
Mx=zeros(n,1);
for i=1:n
    Mx(i,1) = (i-c)*dx;
end
My=zeros(n,1);
for i=1:n
    My(i,1) = (c-i)*dy;
end

r=ro;
for i=1:Nh 
    t = i*dh-dh/2;
    x = ro*cos(t);
    y = ro*sin(t);
    [wx,wxi] = min(abs(Mx-x));
    [wy,wyi] = min(abs(My-y));
    wx = wx+x;
    wy = wy+y;
    if x > wx
        x1 = wx;
        x2 = wx+dx;
        y1 = wy-dy/2;
        y2 = wy+dy/2;
        oyx1y1 = oy(wyi+1,wxi);
        oyx2y1 = oy(wyi+1,wxi+1);
        oyx1y2 = oy(wyi,wxi);
        oyx2y2 = oy(wyi,wxi+1);
        foy = 1/(x2-x1)/(y2-y1)*[x2-x x-x1]*[oyx1y1 oyx1y2; oyx2y1 oyx2y2]*[y2-y; y-y1];
    elseif x == wx
        y1 = wy-dy/2;
        y2 = wy+dy/2;
        oyy1 = oy(wyi+1,wxi);
        oyy2 = oy(wyi,wxi);
        foy = (oyy2-oyy1)/(y2-y1)*(y-y1)+oyy1;
    else
        x1 = wx-dx;
        x2 = wx;
        y1 = wy-dy/2;
        y2 = wy+dy/2;
        oyx1y1 = oy(wyi+1,wxi-1);
        oyx2y1 = oy(wyi+1,wxi);
        oyx1y2 = oy(wyi,wxi-1);
        oyx2y2 = oy(wyi,wxi);
        foy = 1/(x2-x1)/(y2-y1)*[x2-x x-x1]*[oyx1y1 oyx1y2; oyx2y1 oyx2y2]*[y2-y; y-y1];
    end
    if y > wy
        x1 = wx-dx/2;
        x2 = wx+dx/2;
        y1 = wy;
        y2 = wy+dy;
        oxx1y1 = ox(wyi,wxi);
        oxx2y1 = ox(wyi,wxi+1);
        oxx1y2 = ox(wyi-1,wxi);
        oxx2y2 = ox(wyi-1,wxi+1);
        fox = 1/(x2-x1)/(y2-y1)*[x2-x x-x1]*[oxx1y1 oxx1y2; oxx2y1 oxx2y2]*[y2-y; y-y1];
    elseif y == wy
        x1 = wx-dx/2;
        x2 = wx+dx/2;
        oxx1 = ox(wyi,wxi);
        oxx2 = ox(wyi,wxi+1);
        fox = (oxx2-oxx1)/(x2-x1)*(x-x1)+oxx1;
    else
        x1 = wx-dx/2;
        x2 = wx+dx/2;
        y1 = wy-dy;
        y2 = wy;
        oxx1y1 = ox(wyi+1,wxi);
        oxx2y1 = ox(wyi+1,wxi+1);
        oxx1y2 = ox(wyi,wxi);
        oxx2y2 = ox(wyi,wxi+1);
        fox = 1/(x2-x1)/(y2-y1)*[x2-x x-x1]*[oxx1y1 oxx1y2; oxx2y1 oxx2y2]*[y2-y; y-y1];
    end
    l = fox*cos(t)+foy*sin(t);
    or(Nr+1,i) = l;
end