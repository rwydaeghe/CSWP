z = [-50:1:50]
figure(2);
hold on
axis equal

for i = 1:length(z)
        plot(z, z(i)*ones(length(z),1));
        plot(z(i)*ones(length(z),1), z);
end
