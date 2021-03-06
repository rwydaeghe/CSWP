x=[-2:0.05:2]
y=[-2:0.05:2]
q = zeros(81,81)

mov=moviein(200);
figure;
dt=0.001
for it = 1:200 
    t = (it-1)*dt;
    disp([num2str(it) '/' num2str(200)]);
    for i=1:81
        for j = 1:81
        r = x(i)^2+y(81+1-j)^2
        theta = atan2(y(81+1-j),x(i))
        theta = theta.*(theta >= 0) + (theta + 2*pi).*(theta < 0);
        q(j,i) = real(analytic(1,0.1,1,20,r,theta)*exp(-1i*340*1*t));
%         q(j,i) = real(exp(-1i*1*r*cos(theta+pi))*exp(-1i*340*1*t));
        end
    end
    pcolor(q);view(0,90);axis equal;shading interp;caxis([-1 1]);hold on;
    xlim([1 81+1]);ylim([1 81+1]);
    mov(it) = getframe;
end