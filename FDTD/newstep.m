function newstep(Nr,Nh,c,dr,dh,dt,roff)
global or oh p
Z=c;
for i=1:Nh
    for j=1:Nr
        r =  roff+dr/2+(j-1)*dr;
        if i == Nh
            p(j,i) = p(j,i) - c^2*dt/r*(1/dr*((r+dr/2)*or(j+1,i)-(r-dr/2)*or(j,i))+1/dh*(oh(j,1)-oh(j,i)));
        else
            p(j,i) = p(j,i) - c^2*dt/r*(1/dr*((r+dr/2)*or(j+1,i)-(r-dr/2)*or(j,i))+1/dh*(oh(j,i+1)-oh(j,i)));
        end
    end
end


for i=1:Nh
    for j=1:Nr+1
        r=roff+(j-1)*dr;
        if j == 1
            %or(j,i) = 1/(1+Z*dt/dr)*((1-Z*dt/dr)*or(j,i)-2*dt/dr*p(1,i));
        elseif j == Nr+1
            %or(j,i) = 1/(1+Z*dt/dr)*((1-Z*dt/dr)*or(j,i)+2*dt/dr*p(j-1,i));
        else
             or(j,i) = or(j,i) - dt*(1/dr*(p(j,i)-p(j-1,i)));
        end
    end
end

for i=1:Nh
    for j=1:Nr
        r =  roff+dr/2+(j-1)*dr;
        if i == 1
            oh(j,i) = oh(j,i) - dt/r/dh*(p(j,i)-p(j,Nh));
        else
            oh(j,i) = oh(j,i) - dt/r/dh*(p(j,i)-p(j,i-1));
        end
    end
end
        
