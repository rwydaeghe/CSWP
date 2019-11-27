clear all
%test
tic
c=340; %geluidssnelheid - speed of sound (wave speed)
dx=0.04; %ruimtelijke discretisatiestap - spatial discretisation step
dy=dx;
a=1;
nx=round(6*a/dx);
if mod(nx,2)==1
    nx = nx+1; %aantal cellen in x richting - number of cells in x direction
end

ny=nx; %aantal cellen in y richting - number of cells in y direction
nc=round(a/dx)+1;
if mod(nc,2)==1
    nc = nc+1;
end

CFL=1; %Courant getal - Courant number

dt=CFL/(c*sqrt((1/dx^2)+(1/dy^2))); %tijdstap - time step

nt=300/CFL; %aantal tijdstappen in simulatie - number of time steps

x_bron=round(nx/2); y_bron=round(ny/2);



x_recorder=x_bron; y_recorder=y_bron+round(2*a/dx); %plaats recorder 1 - location receiver 1
x_ref=x_bron;y_ref=y_bron+round(2.5*a/dx); %plaats referentie 1 - location reference receiver 1

x_recorder2=x_bron; y_recorder2=y_bron-round(2*a/dx); %plaats recorder 2 - location receiver 2
x_ref2=x_bron;y_ref2=y_bron-round(2.5*a/dx); %plaats referentie 2 - location reference receiver 2

x_recorder3=x_bron-round(2*a/dx); y_recorder3=y_bron; %plaats recorder 3 - location receiver 3
x_ref3=x_bron-round(2.5*a/dx);y_ref3=y_bron; %plaats referentie 3 - location reference receiver 3

x_recorder4=x_bron+round(2*a/dx); y_recorder4=y_bron; %plaats recorder 3 - location receiver 3
x_ref4=x_bron+round(2.5*a/dx);y_ref4=y_bron; %plaats referentie 3 - location reference receiver 3


A=1;k=20;t0=2.5E-2;sigma=5E-5;

global ox oy p
ox = zeros(ny, nx+1); oy = zeros(ny+1, nx); p = zeros(ny, nx);

%film
%movie
mov=moviein(nt);

%initialisatie tijdsreeks recorders
%initialisation time series receivers
recorder = zeros(nt,1);
recorder_ref = zeros(nt,1);

recorder2 = zeros(nt,1);
recorder2_ref = zeros(nt,1);

recorder3 = zeros(nt,1);
recorder3_ref = zeros(nt,1);

recorder4 = zeros(nt,1);
recorder4_ref = zeros(nt,1);

bront = zeros(nt,1);
tijdreeks=zeros(nt,1);
bron=0;

figure;

for it=1:nt
   t = (it-1)*dt; tijdreeks(it)=t;
	disp([num2str(it) '/' num2str(nt)]);
  
   bron=A*sin(2*k*c*(t)); %bron updaten bij nieuw tijd - update source for new time
   
   p(:,1) = p(:,1)+bron; %druk toevoegen bij de drukvergelijking op bronlocatie - adding source term to propagation
   
   step_SIT_SIP_impedance(nx,ny,c,dx,dy,dt,a,nc)   %propagatie over 1 tijdstap - propagate over one time step
  
    recorder(it) = p(x_recorder,y_recorder); %druk opnemen op recorders en referentieplaatsen - store p field at receiver locations
    recorder_ref(it) = p(x_ref,y_ref);
    
   	recorder2(it) = p(x_recorder2,y_recorder2);
    recorder2_ref(it) = p(x_ref2,y_ref2);
    
    recorder3(it) = p(x_recorder3,y_recorder3);
    recorder3_ref(it) = p(x_ref3,y_ref3);
    
    recorder4(it) = p(x_recorder4,y_recorder4);
    recorder4_ref(it) = p(x_ref4,y_ref4);
   
   %voorstellen drukveld
   %presenting the p field   
   pcolor(p);view(0,90);axis equal;shading interp;caxis([-0.02*A 0.02*A]);title([num2str(it) '/' num2str(nt)]);hold on;
   xlim([1 nx+1]);ylim([1 ny+1]);
   plot(x_bron,y_bron,'ks'); plot(x_recorder,y_recorder,'ro');plot(x_ref,y_ref,'ko')
%    plot(x_recorder2,y_recorder2,'ro');plot(x_ref2,y_ref2,'ko');hold off;
   mov(it) = getframe; %wegcommentarieren voor simulatie vlugger te laten lopen - if this line is removed simulation runs much faster
end

n_of_samples=8192;

post_Afout_Pfout(a,dx,dy,c,dt,x_ref,x_recorder,y_ref,y_recorder,x_bron,y_bron,recorder,recorder_ref,n_of_samples,tijdreeks,k)
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref2,x_recorder2,y_ref2,y_recorder2,x_bron,y_bron,recorder2,recorder2_ref,n_of_samples,tijdreeks,k)
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref3,x_recorder3,y_ref3,y_recorder3,x_bron,y_bron,recorder3,recorder3_ref,n_of_samples,tijdreeks,k)
% post_Afout_Pfout(a,dx,dy,c,dt,x_ref4,x_recorder4,y_ref4,y_recorder4,x_bron,y_bron,recorder4,recorder4_ref,n_of_samples,tijdreeks,k)
toc
