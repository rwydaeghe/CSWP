%% Set paramaters

A=1; %amplitude plane wave
k=1; %wave number plane wave
c=340; %speed of sound
Z=c; %surface impedance for non-reflecting boundary
a=1; %radius of circle

dx=0.3; %spatial discretisation step
dy=dx;
nx=round(6*a/dx); %nx numbers of cells in x direction for total field, chosen to be even and large enough to make comparison to analytic solution
if mod(nx,2)==1
    nx = nx+1;
end

ny=nx; %ny chosen equal to nx
npml=30; %number of cells pml
km=10000; %max k for pml
nsf=20; %number of cells sf region
n=nx+2*nsf+2*npml; %total number of cells in x or y direction
CFL=0.1; %Courant number
dt=CFL/(c*sqrt((1/dx^2)+(1/dy^2))); %time step
nt=600/CFL; %number of time steps

%% Set px py ox oy matrices

global ox oy px py 
ox = zeros(n, n+1); oy = zeros(n+1, n); px = zeros(n, n); py = zeros(n, n);
%first index y direction, second index x direction

%% Set reference array induced ox and px (plane wave coming from the left side)

global oxref pxref
oxref = zeros(1,nx+1); pxref = zeros(1,nx);

%% Set circle

global mx my area lxl lxr lyu lyd
positions(a,dx,dy,n);
areas(a,dx,dy,n);
lengths(a,dx,dy,n);

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

%% Set control points

x_bron=round(n/2); y_bron=round(n/2);

x_recorder=x_bron; y_recorder=y_bron+round(2*a/dx); %plaats recorder 1 - location receiver 1
x_ref=x_bron;y_ref=y_bron+round(2.5*a/dx); %plaats referentie 1 - location reference receiver 1

x_recorder2=x_bron; y_recorder2=y_bron-round(2*a/dx); %plaats recorder 2 - location receiver 2
x_ref2=x_bron;y_ref2=y_bron-round(2.5*a/dx); %plaats referentie 2 - location reference receiver 2

x_recorder3=x_bron-round(2*a/dx); y_recorder3=y_bron; %plaats recorder 3 - location receiver 3
x_ref3=x_bron-round(2.5*a/dx);y_ref3=y_bron; %plaats referentie 3 - location reference receiver 3

x_recorder4=x_bron+round(2*a/dx); y_recorder4=y_bron; %plaats recorder 3 - location receiver 3
x_ref4=x_bron+round(2.5*a/dx);y_ref4=y_bron; %plaats referentie 3 - location reference receiver 3


recorder = zeros(nt,1);
recorder_ref = zeros(nt,1);

recorder2 = zeros(nt,1);
recorder2_ref = zeros(nt,1);

recorder3 = zeros(nt,1);
recorder3_ref = zeros(nt,1);

recorder4 = zeros(nt,1);
recorder4_ref = zeros(nt,1);

tijdreeks=zeros(nt,1);

%% Movie

mov=moviein(nt);
figure;

%% Update in time
for it = 1:nt
    t = (it-1)*dt; tijdreeks(it)=t;
    disp([num2str(it) '/' num2str(nt)]);
    
    %% Update incident field
    
    
    pxref(1,1:nx) = pxref(1,1:nx) - c^2*dt/dx*(oxref(1,2:nx+1)-oxref(1,1:nx));
    oxref(1,1) = A*sin(k*c*t);
    oxref(1,2:nx) = oxref(1,2:nx) - dt/dx*(pxref(1,2:nx)-pxref(1,1:nx-1));
    oxref(1,nx+1) = 1/(1+c*dt/dx)*((1-c*dt/dx)*oxref(1,nx+1)+2*dt/dx*pxref(1,nx));
    %     plot(1:dx:1+nx*dx,oxref)
    %    mov(it) = getframe;

    %% Update total p fields
    
    px = ((1-kpx*dt/2).*px - c^2*dt/dx*(ox(:,2:n+1)-ox(:,1:n)))./(1+kpx*dt/2);
    py = ((1-kpy*dt/2).*py - c^2*dt/dy*(oy(1:n,:)-oy(2:n+1,:)))./(1+kpy*dt/2);
%     
    %% Update tf/sf boundary for px (no py component incident field) (px lies in scattered field region)

    px(npml+nsf+1:n-npml-nsf+1,npml+nsf) = px(npml+nsf+1:n-npml-nsf+1,npml+nsf) - c^2*dt/dx*(-oxref(1,1)); %substract (-) oxref from total field on the right (+)
    px(npml+nsf+1:n-npml-nsf+1,n-npml-nsf+1) = px(npml+nsf+1:n-npml-nsf+1,n-npml-nsf+1) - c^2*dt/dx*(oxref(1,nx+1)); %substract (-) oxref from total field on the left(-)
    
    %% Update cells intersecting circle
    lengthm = size(mx,2);
    tol = 0.00000000001;
    for v=1:lengthm
        i = mx(1,v);
        j = my(1,v);
         if area(1,v) < tol
            px(j,i) = 0;
            py(j,i) = 0;
         else
            px(j,i) = px(j,i)+c^2*dt/dx*(ox(j,i+1)-ox(j,i))-c^2*dt/area(1,v)*(lxr(1,v)*ox(j,i+1)-lxl(1,v)*ox(j,i));
            py(j,i) = py(j,i)+c^2*dt/dy*(oy(j,i)-oy(j+1,i))-c^2*dt/area(1,v)*(lyu(1,v)*oy(j,i)-lyd(1,v)*oy(j+1,i));
         end
    end
    %% Update total o fields
    ox(1:n,2:n) = ((1-kox(1:n,2:n)*dt/2).*ox(1:n,2:n) - dt/dx*(px(1:n,2:n)+py(1:n,2:n)-px(1:n,1:n-1)-py(1:n,1:n-1)))./(1+kox(1:n,2:n)*dt/2);
    ox(1:n,1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(1:n,1)-2*dt/dx*(px(1:n,1)+py(1:n,1)));
    ox(1:n,n+1) = 1/(1+Z*dt/dx)*((1-Z*dt/dx)*ox(1:n,n+1)+2*dt/dx*(px(1:n,n)+py(1:n,n)));
    
    oy(2:n,1:n) = ((1-koy(2:n,1:n)*dt/2).*oy(2:n,1:n) - dt/dy*(px(1:n-1,1:n)+py(1:n-1,1:n)-px(2:n,1:n)-py(2:n,1:n)))./(1+koy(2:n,1:n)*dt/2);
      oy(1,1:n) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(1,1:n)+2*dt/dy*(px(1,1:n)+py(1,1:n)));
      oy(n+1,1:n) = 1/(1+Z*dt/dy)*((1-Z*dt/dy)*oy(n+1,1:n)-2*dt/dy*(px(n,1:n)+py(n,1:n)));
%      
    %% Adjust total o fields
    tol = 0.00000000001;
    for v=1:lengthm
        i = mx(1,v);
        j = my(1,v);
        if area(1,v) < tol
            ox(j,i) = 0;
            ox(j,i+1) = 0;
            oy(j,i) = 0;
            oy(j+1,i) = 0;
        end
        if abs(lxl(1,v)) < tol
            ox(j,i) = 0;
        end
        if abs(lxr(1,v)) < tol
            ox(j,i+1) = 0;
        end
        if abs(lyu(1,v)) < tol
            oy(j,i) = 0;
        end
        if abs(lyd(1,v)) < tol
            oy(j+1,i) = 0;
        end
    end
    %% Update tf/sf boundary for ox and oy (both in total field)
    
    ox(npml+nsf+1:n-npml-nsf,npml+nsf+1) = ox(npml+nsf+1:n-npml-nsf,npml+nsf+1) - dt/dx*(-pxref(1,1)); %add (+) incident field to the left (-)
    ox(npml+nsf+1:n-npml-nsf,n-npml-nsf+1) = ox(npml+nsf+1:n-npml-nsf,n-npml-nsf+1) - dt/dx*(+pxref(1,nx)); %add (+) incident field to the right (+)
    
    oy(npml+nsf+1,npml+nsf+1:n-npml-nsf) = oy(npml+nsf+1,npml+nsf+1:n-npml-nsf) - dt/dy*(+pxref(1,:)); %add (+) incident field above (+)
    oy(n-npml-nsf+1,npml+nsf+1:n-npml-nsf) = oy(n-npml-nsf+1,npml+nsf+1:n-npml-nsf) - dt/dy*(-pxref(1,:)); %add (+) incident field below (-)
    
     %% Set control points
    recorder(it) = px(x_recorder,y_recorder)+py(x_recorder,y_recorder); %druk opnemen op recorders en referentieplaatsen - store p field at receiver locations
    recorder_ref(it) = px(x_ref,y_ref)+py(x_ref,y_ref);
    
   	recorder2(it) = px(x_recorder2,y_recorder2)+py(x_recorder2,y_recorder2);
    recorder2_ref(it) = px(x_ref2,y_ref2)+py(x_ref2,y_ref2);
    
    recorder3(it) = px(x_recorder3,y_recorder3)+py(x_recorder3,y_recorder3);
    recorder3_ref(it) = px(x_ref3,y_ref3)+py(x_ref3,y_ref3);
    
    recorder4(it) = px(x_recorder4,y_recorder4)+py(x_recorder4,y_recorder4);
    recorder4_ref(it) = px(x_ref4,y_ref4)+py(x_ref4,y_ref4);
    %% Presenting the p field
%     
%        pcolor(px+py);view(0,90);axis equal;shading interp;caxis([-340*A 340*A]);title([num2str(it) '/' num2str(nt)]);hold on;
%        xlim([1 n+1]);ylim([1 n+1]);
%        mov(it) = getframe; %if this line is removed simulation runs much faster
%     
%     
end

pcolor(px+py);view(0,90);axis equal;shading interp;caxis([-340*A 340*A]);title([num2str(it) '/' num2str(nt)]);hold on;
xlim([1 n+1]);ylim([1 n+1]);
n_of_samples=8192;
% % 
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref,x_recorder,y_ref,y_recorder,x_bron,y_bron,recorder,recorder_ref,n_of_samples,tijdreeks,k)
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref2,x_recorder2,y_ref2,y_recorder2,x_bron,y_bron,recorder2,recorder2_ref,n_of_samples,tijdreeks,k)
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref3,x_recorder3,y_ref3,y_recorder3,x_bron,y_bron,recorder3,recorder3_ref,n_of_samples,tijdreeks,k)
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref4,x_recorder4,y_ref4,y_recorder4,x_bron,y_bron,recorder4,recorder4_ref,n_of_samples,tijdreeks,k)

