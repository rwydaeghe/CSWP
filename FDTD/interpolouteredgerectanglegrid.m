function interpolouteredgerectanglegrid(nr, nh, n, dx, dy, dh, dr, ro)
global q
global ox oy
global or oh
dr=0.1;
ro=1;
nr=5;
Mr = [ro+dr/2:dr:ro+nr*dr];
Mh = [0+dh/2:dh:2*pi];
c = n/2+0.5;
for j=1:n
    for i=1:n-1
        if q(j,i) ~= q(j,i+1)
            x = (i-c)*dx+dx/2;
            y = (c-j)*dy;
            r = sqrt(x^2+y^2);
            h = atan2(y, x); 
            h = h.*(h >= 0) + (h + 2*pi).*(h < 0);
            [wr,wri] = min(abs(Mr-r));
            [wh,whi] = min(abs(Mh-h));
            if h < wh
                h1 = wh;
                h2 = wh-dh;
                r1 = wr-dr/2;
                r2 = wr+dr/2;
                orh1r1 = or(wri,whi);
                orh1r2 = or(wri+1,whi);
                if whi == 1
                    orh2r1 = or(wri,nh);
                    orh2r2 = or(wri+1,nh);
                else
                    orh2r1 = or(wri,whi-1);
                    orh2r2 = or(wri+1,whi-1);
                end
                ior = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[orh1r1 orh1r2; orh2r1 orh2r2]*[r2-r; r-r1];      
            else
                h1 = wh+dh;
                h2 = wh;
                r1 = wr-dr/2;
                r2 = wr+dr/2;
                orh2r1 = or(wri,whi);
                orh2r2 = or(wri+1,whi);
                if whi == nh
                    orh1r1 = or(wri,1);
                    orh1r2 = or(wri+1,1);
                else
                    orh1r1 = or(wri,whi+1);
                    orh1r2 = or(wri+1,whi+1);
                end
                ior = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[orh1r1 orh1r2; orh2r1 orh2r2]*[r2-r; r-r1];
            end
            if r > wr
                h1 = wh+dh/2;
                h2 = wh-dh/2;
                r1 = wr;
                r2 = wr+dr;
                ohh2r1 = oh(wri,whi);
                ohh2r2 = oh(wri+1,whi);
                if whi == nh
                    ohh1r1 = oh(wri,1);
                    ohh1r2 = oh(wri+1,1);
                else
                    ohh1r1 = oh(wri,whi+1);
                    ohh1r2 = oh(wri+1,whi+1);
                end
                ioh = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[ohh1r1 ohh1r2; ohh2r1 ohh2r2]*[r2-r; r-r1];
            else
                h1 = wh+dh/2;
                h2 = wh-dh/2;
                r1 = wr-dr;
                r2 = wr;
                ohh2r1 = oh(wri-1,whi);
                ohh2r2 = oh(wri,whi);
                if whi == nh
                    ohh1r1 = oh(wri-1,1);
                    ohh1r2 = oh(wri,1);
                else
                    ohh1r1 = oh(wri-1,whi+1);
                    ohh1r2 = oh(wri,whi+1);
                end
                ioh = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[ohh1r1 ohh1r2; ohh2r1 ohh2r2]*[r2-r; r-r1];
            end
            ox(j,i+1) = ior*cos(h)-ioh*sin(h);
        end
    end
end
for j=1:n-1
    for i=1:n
        if q(j,i) ~= q(j+1,i)
            x = (i-c)*dx;
            y = (c-j)*dy+dy/2;
            r = sqrt(x^2+y^2);
            h = atan2(y, x); 
            h = h.*(h >= 0) + (h + 2*pi).*(h < 0);
            [wr,wri] = min(abs(Mr-r));
            [wh,whi] = min(abs(Mh-h));
            if h < wh
                h1 = wh;
                h2 = wh-dh;
                r1 = wr-dr/2;
                r2 = wr+dr/2;
                orh1r1 = or(wri,whi);
                orh1r2 = or(wri+1,whi);
                if whi == 1
                    orh2r1 = or(wri,nh);
                    orh2r2 = or(wri+1,nh);
                else
                    orh2r1 = or(wri,whi-1);
                    orh2r2 = or(wri+1,whi-1);
                end
                ior = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[orh1r1 orh1r2; orh2r1 orh2r2]*[r2-r; r-r1];      
            else
                h1 = wh+dh;
                h2 = wh;
                r1 = wr-dr/2;
                r2 = wr+dr/2;
                orh2r1 = or(wri,whi);
                orh2r2 = or(wri+1,whi);
                if whi == nh
                    orh1r1 = or(wri,1);
                    orh1r2 = or(wri+1,1);
                else
                    orh1r1 = or(wri,whi+1);
                    orh1r2 = or(wri+1,whi+1);
                end
                ior = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[orh1r1 orh1r2; orh2r1 orh2r2]*[r2-r; r-r1];
            end
            if r > wr
                h1 = wh+dh/2;
                h2 = wh-dh/2;
                r1 = wr;
                r2 = wr+dr;
                ohh2r1 = oh(wri,whi);
                ohh2r2 = oh(wri+1,whi);
                if whi == nh
                    ohh1r1 = oh(wri,1);
                    ohh1r2 = oh(wri+1,1);
                else
                    ohh1r1 = oh(wri,whi+1);
                    ohh1r2 = oh(wri+1,whi+1);
                end
                ioh = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[ohh1r1 ohh1r2; ohh2r1 ohh2r2]*[r2-r; r-r1];
            else
                h1 = wh+dh/2;
                h2 = wh-dh/2;
                r1 = wr-dr;
                r2 = wr;
                ohh2r1 = oh(wri-1,whi);
                ohh2r2 = oh(wri,whi);
                if whi == nh
                    ohh1r1 = oh(wri-1,1);
                    ohh1r2 = oh(wri,1);
                else
                    ohh1r1 = oh(wri-1,whi+1);
                    ohh1r2 = oh(wri,whi+1);
                end
                ioh = 1/(h2-h1)/(r2-r1)*[h2-h h-h1]*[ohh1r1 ohh1r2; ohh2r1 ohh2r2]*[r2-r; r-r1];
            end
            oy(j,i+1) = ior*sin(h)+ioh*cos(h);
        end
    end
end
    
