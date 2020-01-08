global or oh p
%% Set paramaters

A=1; %amplitude plane wave
k=5; %wave number plane wave
c=340; %speed of sound
Z=c; %surface impedance for non-reflecting boundary
a=1; %radius of circle

dx=0.1; %spatial discretisation step
dy=dx;
nx=round(6*a/dx); %nx numbers of cells in x direction for total field, chosen to be even and large enough to make comparison to analytic solution
if mod(nx,2)==1
    nx = nx+1;
end
ny=nx; %ny chosen equal to nx

dr=dx; %spatial discretisation step in r direction in polar grid (dr = dx, cells need to have approx the same size)
nr=20; %number of cells in r direction (only 4)
roff=a; %inner edge polar grid (starts at the edge of the circle)
ro=roff+nr*dr; %outer edge polar grid
dh=dr/roff; %spatial discretisation step in theta direction
nh=2*round(pi/dh); %number of cells in theta direction
dh=2*pi/nh; %dh recalculated so the whole 0..2pi interval is covered

npml=30; %number of cells pml
km=1000; %max k for pml
nsf=20; %number of cells sf region
n=nx+2*nsf+2*npml; %total number of cells in x or y direction
CFL=1; %Courant number
dt=CFL/(c*sqrt((1/dx^2)+(1/dy^2))); %time step
nt=200/CFL; %number of time steps

R=[roff+dr/2:dr:ro-dr/2];
T=[dh/2:dh:2*pi-dh/2]*180/pi;
for it = 1:nt
      t = (it-1)*dt;
    disp([num2str(it) '/' num2str(nt)]);
    newstep(nr,nh,c,dr,dh,dt,roff)
     polarPcolor(R,T,p); shading interp; caxis([-340*A 340*A]); title([num2str(it) '/' num2str(nt)]);
    %pcolor(px+py);view(0,90);axis equal;shading interp;caxis([-340*A 340*A]);title([num2str(it) '/' num2str(nt)]);hold on;
    %xlim([1 n+1]);ylim([1 n+1]);
    mov(it) = getframe; %if this line is removed simulation runs much faster
end