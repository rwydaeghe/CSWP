function lengths(a,dx,dy,n)
cc = n/2+0.5;
global mx my lxl lxr lyu lyd
lxl=[];
lxr=[];
lyu=[];
lyd=[];
for v=1:size(mx,2)
    i = mx(1,v);
    j = my(1,v);
    x = (i-cc)*dx;
    y = (cc-j)*dy;
    x1=x-dx/2;
    x2=x+dx/2;
    y1=y-dy/2;
    y2=y+dy/2;
    sl = dx/1000000;
    ttl = dx/sl+1;
    l = [x1:sl:x2];
    tll = 0.1;
    tlr = 0.1;
    tlu = 0.1;
    tld = 0.1;
    ytl=sqrt(a^2-x1^2);
    if x1^2+y1^2<=a^2 && x1^2+y2^2<=a^2
        tll = 0;
    elseif x1^2+y1^2<=a^2 || x1^2+y2^2<=a^2
        if y1 <= ytl && ytl <= y2
            tll = y2-ytl;
        elseif y1 <= -ytl && -ytl <= y2
            tll = -ytl-y1;
        end
    end
    
    ytr=sqrt(a^2-x2^2);
    if x2^2+y1^2<=a^2 && x2^2+y2^2<=a^2
        tlr = 0;
    elseif x2^2+y1^2<=a^2 || x2^2+y2^2<=a^2
        if y1 <= ytr && ytr <= y2
            tlr = y2-ytr;
        elseif y1 <= -ytr && -ytr <= y2
            tlr = -ytr-y1;
        end
    end
    
    xtu=sqrt(a^2-y2^2);
    if x1^2+y2^2<=a^2 && x2^2+y2^2<=a^2
        tlu = 0;
    elseif x1^2+y2^2<=a^2 || x2^2+y2^2<=a^2
        if x1 <= xtu && xtu <= x2
            tlu = x2-xtu;
        elseif x1 <= -xtu && -xtu <= x2
            tlu = -xtu-x1;
        end
    end
    
    xtd=sqrt(a^2-y1^2);
    if x1^2+y1^2<=a^2 && x2^2+y1^2<=a^2
        tld = 0;
    elseif x1^2+y1^2<=a^2 || x2^2+y1^2<=a^2
        if x1 <= xtd && xtd <= x2
            tld = x2-xtd;
        elseif x1 <= -xtd && -xtd <= x2
            tld = -xtd-x1;
        end
    end
    lxl = [lxl tll];
    lxr = [lxr tlr];
    lyu = [lyu tlu];
    lyd = [lyd tld];
end

