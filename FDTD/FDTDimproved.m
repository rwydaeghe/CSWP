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
nr=6; %number of cells in r direction (only 4)
roff=a; %inner edge polar grid (starts at the edge of the circle)
ro=roff+nr*dr; %outer edge polar grid
dh=dr/roff; %spatial discretisation step in theta direction
nh=round(2*pi/dh); %number of cells in theta direction
dh=2*pi/nh; %dh recalculated so the whole 0..2pi interval is covered

npml=30; %number of cells pml
km=1000; %max k for pml
nsf=20; %number of cells sf region
n=nx+2*nsf+2*npml; %total number of cells in x or y direction
CFL=1; %Courant number
dt=CFL/(c*sqrt((1/dx^2)+(1/dy^2))); %time step
nt=400/CFL; %number of time steps

%% Set px py ox oy matrices

global ox oy px py 
ox = zeros(n, n+1); oy = zeros(n+1, n); px = zeros(n, n); py = zeros(n, n);
%first index y direction, second index x direction

%% Set p or oh matrices

global p or oh 
or = zeros(nr+1, nh); oh = zeros(nr, nh); p = zeros(nr, nh);
%first index y direction, second index x direction

%% Set reference array induced ox and px (plane wave coming from the left side)

global oxref pxref
oxref = zeros(1,nx+1); pxref = zeros(1,nx);

%% Set circle

global q qx qy %matrices which after elementwise multiplication with the px py ox and oy matrix induce the boundary condition for the fields, using the staircase method to get the circle geometry
q = ones(n,n);
qx = ones(n,n+1);
qy = ones(n+1,n);
circlewithpolar(dx,dy,a,n);

%% Set k values for the whole region (to have pml)

kpx = zeros(1,n);
kox = zeros(1,n+1);

for i = 1:n
    if i <= npml
        kpx(i) = km*((npml-i+1/2)/npml)^4;
    elseif i >= n-npml+1
        kpx(i) = km*((i-n+npml-1/2)/npml)^4;
    else
        kpx(i) = 0;
    end
end
for i = 1:n+1
    if i <= npml
        kox(i) = km*((npml-i)/npml)^4;
    elseif i >= n-npml+1
        kox(i) = km*((i-n+npml)/npml)^4;
    else
        kox(i) = 0;
    end
end

kpx = repmat(kpx,n,1);
kpy = kpx';
kox = repmat(kox,n,1);
koy = kox';


%% Movie

mov=moviein(nt);
figure;

%% Update in time
for it = 1:nt
    t = (it-1)*dt; %tijdreeks(it)=t;
    disp([num2str(it) '/' num2str(nt)]);
    
    %% Update incident field
    
    
    pxref(1,1:nx) = pxref(1,1:nx) - c^2*dt/dx*(oxref(1,2:nx+1)-oxref(1,1:nx));
    oxref(1,1) = A*sin(k*c*t);
    oxref(1,2:nx) = oxref(1,2:nx) - dt/dx*(pxref(1,2:nx)-pxref(1,1:nx-1));
    oxref(1,nx+1) = 1/(1+c*dt/dx)*((1-c*dt/dx)*oxref(1,nx+1)+2*dt/dx*pxref(1,nx));
    %     plot(1:dx:1+nx*dx,oxref)
    %    mov(it) = getframe;
    
    %% Interpolate for boundary conditions on rectangular grid
    interpolouteredgerectanglegrid(nr,nh,n,dx,dy,dh,dr,ro)
    

    %% Update total p fields rectangular grid

    px = ((1-kpx*dt/2).*px - c^2*dt/dx*(ox(:,2:n+1)-ox(:,1:n)))./(1+kpx*dt/2);
    py = ((1-kpy*dt/2).*py - c^2*dt/dy*(oy(2:n+1,:)-oy(1:n,:)))./(1+kpy*dt/2);
    
    %% Update tf/sf boundary for px (no py component incident field) (px lies in scattered field region)

    px(npml+nsf+1:n-npml-nsf+1,npml+nsf) = px(npml+nsf+1:n-npml-nsf+1,npml+nsf) - c^2*dt/dx*(-oxref(1,1)); %substract (-) oxref from total field on the right (+)
    px(npml+nsf+1:n-npml-nsf+1,n-npml-nsf+1) = px(npml+nsf+1:n-npml-nsf+1,n-npml-nsf+1) - c^2*dt/dx*(oxref(1,nx+1)); %substract (-) oxref from total field on the left(-)

    
    %% Update total o fields rectangular grid
    
    ox(1:n,2:n) = ((1-kox(1:n,2:n)*dt/2).*ox(1:n,2:n) - dt/dx*(px(1:n,2:n)+py(1:n,2:n)-px(1:n,1:n-1)-py(1:n,1:n-1)))./(1+kox(1:n,2:n)*dt/2);
    ox(1:n,1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(1:n,1)-2*dt/dx*(px(1:n,1)+py(1:n,1)));
    ox(1:n,n+1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(1:n,n+1)+2*dt/dx*(px(1:n,n)+py(1:n,n)));
    
     oy(2:n,1:n) = ((1-koy(2:n,1:n)*dt/2).*oy(2:n,1:n) - dt/dy*(px(2:n,1:n)+py(2:n,1:n)-px(1:n-1,1:n)-py(1:n-1,1:n)))./(1+koy(2:n,1:n)*dt/2);
     oy(1,1:n) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(1,1:n)-2*dt/dy*(px(1,1:n)+py(1,1:n)));
     oy(n+1,1:n) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(n+1,1:n)+2*dt/dy*(px(n,1:n)+py(n,1:n)));
     
    
    %% Update tf/sf boundary for ox and oy (both in total field)
    
    ox(npml+nsf+1:n-npml-nsf,npml+nsf+1) = ox(npml+nsf+1:n-npml-nsf,npml+nsf+1) - dt/dx*(-pxref(1,1)); %add (+) incident field to the left (-)
    ox(npml+nsf+1:n-npml-nsf,n-npml-nsf+1) = ox(npml+nsf+1:n-npml-nsf,n-npml-nsf+1) - dt/dx*(+pxref(1,nx)); %add (+) incident field to the right (+)
    
    oy(npml+nsf+1,npml+nsf+1:n-npml-nsf) = oy(npml+nsf+1,npml+nsf+1:n-npml-nsf) - dt/dy*(-pxref(1,:)); %add (+) incident field below (-)
    oy(n-npml-nsf+1,npml+nsf+1:n-npml-nsf) = oy(n-npml-nsf+1,npml+nsf+1:n-npml-nsf) - dt/dy*(pxref(1,:)); %add (+) incident field above (+)
    
    %% Ensure boundary cylinder
    
     %ox(npml+nsf+1:npml+nsf+ny,npml+nsf+1:npml+nsf+nx+1) = ox(npml+nsf+1:npml+nsf+ny,npml+nsf+1:npml+nsf+nx+1).*qx;
     %oy(npml+nsf+1:npml+nsf+ny+1,npml+nsf+1:npml+nsf+nx) = oy(npml+nsf+1:npml+nsf+ny+1,npml+nsf+1:npml+nsf+nx).*qy;
     %px(npml+nsf+1:npml+nsf+ny,npml+nsf+1:npml+nsf+nx) = px(npml+nsf+1:npml+nsf+ny,npml+nsf+1:npml+nsf+nx).*q;
     %py(npml+nsf+1:npml+nsf+ny,npml+nsf+1:npml+nsf+nx) = py(npml+nsf+1:npml+nsf+ny,npml+nsf+1:npml+nsf+nx).*q;
     
     ox = ox.*qx;
     oy = oy.*qy;
     px = px.*q;
     py = py.*q;
     %% Interpolate for boundary conditions on polar grid (outer edge)
     interpolouteredgecircle(nr,nh,ro,dr,dh,n,dx,dy)
     
     %% FDTD on polar grid
     newstep(nr,nh,c,dr,dh,dt,roff)
     
    %% Presenting the p field
    
    pcolor(px+py);view(0,90);axis equal;shading interp;caxis([-340*A 340*A]);title([num2str(it) '/' num2str(nt)]);hold on;
    xlim([1 n+1]);ylim([1 n+1]);
    mov(it) = getframe; %if this line is removed simulation runs much faster
    
    
end
    







