function positions16(a,dx,dy,n)
tol = 0.000000000000001;
cc = n/2+0.5;
global mx16 my16 locx locy
mx16=[];
my16=[];
locx = [];
locy = [];
cx1 = -9*a/2;
cx2 = -3*a/2;
cx3 = 3*a/2;
cx4 = 9*a/2;
cy1 = -9*a/2;
cy2 = -3*a/2;
cy3 = 3*a/2;
cy4 = 9*a/2;

for i = 1:n
    for j = 1:n
        x = (i-cc)*dx;
        y = (cc-j)*dy;
        x1=x-dx/2;
        x2=x+dx/2;
        y1=y-dy/2;
        y2=y+dy/2;
        %case1
        if (x1-cx1)^2+(y1-cy1)^2 <= a^2 || (x2-cx1)^2+(y1-cy1)^2 <= a^2 || (x1-cx1)^2+(y2-cy1)^2 <= a^2 || (x2-cx1)^2+(y2-cy1)^2 <= a^2
            if (x1-cx1)^2+(y1-cy1)^2 >= a^2 & (x2-cx1)^2+(y1-cy1)^2 >= a^2 & (x1-cx1)^2+(y2-cy1)^2 >= a^2 & (x2-cx1)^2+(y2-cy1)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx1];
                locy = [locy cy1];
            end
        end
        
        %case2
        if (x1-cx1)^2+(y1-cy2)^2 <= a^2 || (x2-cx1)^2+(y1-cy2)^2 <= a^2 || (x1-cx1)^2+(y2-cy2)^2 <= a^2 || (x2-cx1)^2+(y2-cy2)^2 <= a^2
            if (x1-cx1)^2+(y1-cy2)^2 >= a^2 & (x2-cx1)^2+(y1-cy2)^2 >= a^2 & (x1-cx1)^2+(y2-cy2)^2 >= a^2 & (x2-cx1)^2+(y2-cy2)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx1];
                locy = [locy cy2];
            end
        end
        
        %case3
        if (x1-cx1)^2+(y1-cy3)^2 <= a^2 || (x2-cx1)^2+(y1-cy3)^2 <= a^2 || (x1-cx1)^2+(y2-cy3)^2 <= a^2 || (x2-cx1)^2+(y2-cy3)^2 <= a^2
            if (x1-cx1)^2+(y1-cy3)^2 >= a^2 & (x2-cx1)^2+(y1-cy3)^2 >= a^2 & (x1-cx1)^2+(y2-cy3)^2 >= a^2 & (x2-cx1)^2+(y2-cy3)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx1];
                locy = [locy cy3];
            end
        end
        
        %case4
        if (x1-cx1)^2+(y1-cy4)^2 <= a^2 || (x2-cx1)^2+(y1-cy4)^2 <= a^2 || (x1-cx1)^2+(y2-cy4)^2 <= a^2 || (x2-cx1)^2+(y2-cy4)^2 <= a^2
            if (x1-cx1)^2+(y1-cy4)^2 >= a^2 & (x2-cx1)^2+(y1-cy4)^2 >= a^2 & (x1-cx1)^2+(y2-cy4)^2 >= a^2 & (x2-cx1)^2+(y2-cy4)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx1];
                locy = [locy cy4];
            end
        end
        
        %case5
        if (x1-cx2)^2+(y1-cy1)^2 <= a^2 || (x2-cx2)^2+(y1-cy1)^2 <= a^2 || (x1-cx2)^2+(y2-cy1)^2 <= a^2 || (x2-cx2)^2+(y2-cy1)^2 <= a^2
            if (x1-cx2)^2+(y1-cy1)^2 >= a^2 & (x2-cx2)^2+(y1-cy1)^2 >= a^2 & (x1-cx2)^2+(y2-cy1)^2 >= a^2 & (x2-cx2)^2+(y2-cy1)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx2];
                locy = [locy cy1];
            end
        end
        
        %case6
        if (x1-cx2)^2+(y1-cy2)^2 <= a^2 || (x2-cx2)^2+(y1-cy2)^2 <= a^2 || (x1-cx2)^2+(y2-cy2)^2 <= a^2 || (x2-cx2)^2+(y2-cy2)^2 <= a^2
            if (x1-cx2)^2+(y1-cy2)^2 >= a^2 & (x2-cx2)^2+(y1-cy2)^2 >= a^2 & (x1-cx2)^2+(y2-cy2)^2 >= a^2 & (x2-cx2)^2+(y2-cy2)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx2];
                locy = [locy cy2];
            end
        end
        
        %case7
        if (x1-cx2)^2+(y1-cy3)^2 <= a^2 || (x2-cx2)^2+(y1-cy3)^2 <= a^2 || (x1-cx2)^2+(y2-cy3)^2 <= a^2 || (x2-cx2)^2+(y2-cy3)^2 <= a^2
            if (x1-cx2)^2+(y1-cy3)^2 >= a^2 & (x2-cx2)^2+(y1-cy3)^2 >= a^2 & (x1-cx2)^2+(y2-cy3)^2 >= a^2 & (x2-cx2)^2+(y2-cy3)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx2];
                locy = [locy cy3];
            end
        end
        
        %case8
        if (x1-cx2)^2+(y1-cy4)^2 <= a^2 || (x2-cx2)^2+(y1-cy4)^2 <= a^2 || (x1-cx2)^2+(y2-cy4)^2 <= a^2 || (x2-cx2)^2+(y2-cy4)^2 <= a^2
            if (x1-cx2)^2+(y1-cy4)^2 >= a^2 & (x2-cx2)^2+(y1-cy4)^2 >= a^2 & (x1-cx2)^2+(y2-cy4)^2 >= a^2 & (x2-cx2)^2+(y2-cy4)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx2];
                locy = [locy cy4];
            end
        end
        
        %case9
        if (x1-cx3)^2+(y1-cy1)^2 <= a^2 || (x2-cx3)^2+(y1-cy1)^2 <= a^2 || (x1-cx3)^2+(y2-cy1)^2 <= a^2 || (x2-cx3)^2+(y2-cy1)^2 <= a^2
            if (x1-cx3)^2+(y1-cy1)^2 >= a^2 & (x2-cx3)^2+(y1-cy1)^2 >= a^2 & (x1-cx3)^2+(y2-cy1)^2 >= a^2 & (x2-cx3)^2+(y2-cy1)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx3];
                locy = [locy cy1];
            end
        end
        
        %case10
        if (x1-cx3)^2+(y1-cy2)^2 <= a^2 || (x2-cx3)^2+(y1-cy2)^2 <= a^2 || (x1-cx3)^2+(y2-cy2)^2 <= a^2 || (x2-cx3)^2+(y2-cy2)^2 <= a^2
            if (x1-cx3)^2+(y1-cy2)^2 >= a^2 & (x2-cx3)^2+(y1-cy2)^2 >= a^2 & (x1-cx3)^2+(y2-cy2)^2 >= a^2 & (x2-cx3)^2+(y2-cy2)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx3];
                locy = [locy cy2];
            end
        end
        
        %case11
        if (x1-cx3)^2+(y1-cy3)^2 <= a^2 || (x2-cx3)^2+(y1-cy3)^2 <= a^2 || (x1-cx3)^2+(y2-cy3)^2 <= a^2 || (x2-cx3)^2+(y2-cy3)^2 <= a^2
            if (x1-cx3)^2+(y1-cy3)^2 >= a^2 & (x2-cx3)^2+(y1-cy3)^2 >= a^2 & (x1-cx3)^2+(y2-cy3)^2 >= a^2 & (x2-cx3)^2+(y2-cy3)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx3];
                locy = [locy cy3];
            end
        end
        
        %case12
        if (x1-cx3)^2+(y1-cy4)^2 <= a^2 || (x2-cx3)^2+(y1-cy4)^2 <= a^2 || (x1-cx3)^2+(y2-cy4)^2 <= a^2 || (x2-cx3)^2+(y2-cy4)^2 <= a^2
            if (x1-cx3)^2+(y1-cy4)^2 >= a^2 & (x2-cx3)^2+(y1-cy4)^2 >= a^2 & (x1-cx3)^2+(y2-cy4)^2 >= a^2 & (x2-cx3)^2+(y2-cy4)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx3];
                locy = [locy cy4];
            end
        end
        
        %case13
        if (x1-cx4)^2+(y1-cy1)^2 <= a^2 || (x2-cx4)^2+(y1-cy1)^2 <= a^2 || (x1-cx4)^2+(y2-cy1)^2 <= a^2 || (x2-cx4)^2+(y2-cy1)^2 <= a^2
            if (x1-cx4)^2+(y1-cy1)^2 >= a^2 & (x2-cx4)^2+(y1-cy1)^2 >= a^2 & (x1-cx4)^2+(y2-cy1)^2 >= a^2 & (x2-cx4)^2+(y2-cy1)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx4];
                locy = [locy cy1];
            end
        end
        
        %case14
        if (x1-cx4)^2+(y1-cy2)^2 <= a^2 || (x2-cx4)^2+(y1-cy2)^2 <= a^2 || (x1-cx4)^2+(y2-cy2)^2 <= a^2 || (x2-cx4)^2+(y2-cy2)^2 <= a^2
            if (x1-cx4)^2+(y1-cy2)^2 >= a^2 & (x2-cx4)^2+(y1-cy2)^2 >= a^2 & (x1-cx4)^2+(y2-cy2)^2 >= a^2 & (x2-cx4)^2+(y2-cy2)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx4];
                locy = [locy cy2];
            end
        end
        
        %case15
        if (x1-cx4)^2+(y1-cy3)^2 <= a^2 || (x2-cx4)^2+(y1-cy3)^2 <= a^2 || (x1-cx4)^2+(y2-cy3)^2 <= a^2 || (x2-cx4)^2+(y2-cy3)^2 <= a^2
            if (x1-cx4)^2+(y1-cy3)^2 >= a^2 & (x2-cx4)^2+(y1-cy3)^2 >= a^2 & (x1-cx4)^2+(y2-cy3)^2 >= a^2 & (x2-cx4)^2+(y2-cy3)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx4];
                locy = [locy cy3];
            end
        end
        
        %case16
        if (x1-cx4)^2+(y1-cy4)^2 <= a^2 || (x2-cx4)^2+(y1-cy4)^2 <= a^2 || (x1-cx4)^2+(y2-cy4)^2 <= a^2 || (x2-cx4)^2+(y2-cy4)^2 <= a^2
            if (x1-cx4)^2+(y1-cy4)^2 >= a^2 & (x2-cx4)^2+(y1-cy4)^2 >= a^2 & (x1-cx4)^2+(y2-cy4)^2 >= a^2 & (x2-cx4)^2+(y2-cy4)^2 >= a^2
            else
                mx16 = [mx16 i];
                my16 = [my16 j];
                locx = [locx cx4];
                locy = [locy cy4];
            end
        end
    end
end



                    
            