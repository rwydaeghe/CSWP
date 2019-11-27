ri=40;
ro=50;
N=10;
dr=(ro-ri)/N;
A=5;
da=2*pi/A;
r=[ri:dr:ro];
a=[0:da:2*pi];
q=[0:0.01:2*pi]
z = [-50:1:50]


figure (1);
hold on
axis equal
for i = 1:length(r)
    plot(r(i)*cos(q),r(i)*sin(q));
end

for i = 1:length(a)
    plot([ri ro]*cos(a(i)),[ri ro]*sin(a(i)));
end
for i = 1:length(z)
        plot(z, z(i)*ones(length(z),1));
        plot(z(i)*ones(length(z),1), z);
end


    



