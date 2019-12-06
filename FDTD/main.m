clear all
directory='c:\oefening' %INSTELLEN ! ! - MODIFY
addpath(directory)

%INITIALISATIE 2D-GRID EN SIMULATIEPARAMETERS-----------------------------
%INITIALISATION 2D-GRID AND SIMULATION PARAMETERS-------------------------

c=340; %geluidssnelheid - speed of sound (wave speed)
dx=0.2; %ruimtelijke discretisatiestap - spatial discretisation step
dy=dx;

nx=250; %aantal cellen in x richting - number of cells in x direction
ny=250; %aantal cellen in y richting - number of cells in y direction

CFL=1; %Courant getal - Courant number

dt=CFL/(c*sqrt((1/dx^2)+(1/dy^2))); %tijdstap - time step

nt=400/CFL; %aantal tijdstappen in simulatie - number of time steps

%locatie bron (midden van grid) en ontvangers
%location of source(central) and receivers
x_bron=round(nx/2); y_bron=round(ny/2);

x_recorder=x_bron+30; y_recorder=y_bron+30; %plaats recorder 1 - location receiver 1
x_ref=x_bron+20;y_ref=y_bron+20; %plaats referentie 1 - location reference receiver 1

x_recorder2=x_bron+30; y_recorder2=y_bron; %plaats recorder 2 - location receiver 2
x_ref2=x_bron+20;y_ref2=y_bron; %plaats referentie 2 - location reference receiver 2

%pulse gegevens 
%source pulse information
A=1;fc=100;t0=2.5E-2;sigma=5E-5;

%initialisatie snelheids- en drukvelden
%initialisation of o and p fields 

n=40;
global ox oy p
ox = zeros(ny, nx+1); oy = zeros(ny+1, nx); p = zeros(ny, nx); 

global oxl oyl plo plp
oxl = zeros(ny, n); oyl = zeros(ny+1, n); plo = zeros(ny, n); plp = zeros(ny, n);

global oxr oyr pro prp
oxr = zeros(ny, n); oyr = zeros(ny+1, n); pro = zeros(ny, n); prp = zeros(ny, n);

global oxu oyu puo pup
oxu = zeros(n, nx+1); oyu = zeros(n, nx); puo = zeros(n, nx); pup = zeros(n, nx);

global oxd oyd pdo pdp
oxd = zeros(n, nx+1); oyd = zeros(n, nx); pdo = zeros(n, nx); pdp = zeros(n, nx);

%film
%movie
mov=moviein(nt);

%initialisatie tijdsreeks recorders
%initialisation time series receivers
recorder = zeros(nt,1);
recorder_ref = zeros(nt,1);

recorder2 = zeros(nt,1);
recorder2_ref = zeros(nt,1);

bront = zeros(nt,1);
tijdreeks=zeros(nt,1);
bron=0;

figure;
a=1;
nc=1;

%TIJDSITERATIE------------------------------------------------------
%TIME ITTERATION----------------------------------------------------
for it=1:nt
   t = (it-1)*dt; tijdreeks(it)=t;
	disp([num2str(it) '/' num2str(nt)]);
  
   bron=A*sin(2*pi*fc*(t-t0))*exp(-((t-t0)^2)/(sigma)); %bron updaten bij nieuw tijd - update source for new time
   
   p(x_bron,y_bron) = p(x_bron,y_bron)+bron; %druk toevoegen bij de drukvergelijking op bronlocatie - adding source term to propagation
   
   step_SIT_SIP_impedance(nx,ny,c,dx,dy,dt,a,nc)   %propagatie over 1 tijdstap - propagate over one time step
  
    recorder(it) = p(x_recorder,y_recorder); %druk opnemen op recorders en referentieplaatsen - store p field at receiver locations
   recorder_ref(it) = p(x_ref,y_ref);
   
  	recorder2(it) = p(x_recorder2,y_recorder2);
   recorder2_ref(it) = p(x_ref2,y_ref2);
   
   %voorstellen drukveld
   %presenting the p field   
   pcolor(pup);view(0,90);axis equal;shading interp;caxis([-0.02*A 0.02*A]);title([num2str(it) '/' num2str(nt)]);hold on;
   xlim([1 nx+1]);ylim([1 ny+1]);
   plot(x_bron,y_bron,'ks');plot(x_recorder,y_recorder,'ro');plot(x_ref,y_ref,'ko')
   plot(x_recorder2,y_recorder2,'ro');plot(x_ref2,y_ref2,'ko');hold off;
   mov(it) = getframe; %wegcommentarieren voor simulatie vlugger te laten lopen - if this line is removed simulation runs much faster
end

%movie(mov); %laten afspelen opgenomen simulatie - play back of stored movie

%NAVERWERKING : BEREKENING FASEFOUT en AMPLITUDEFOUT---------------------------------
%POST PROCESSING : CALCULATE PHASE and AMPLITUDE ERROR-------------------------------
n_of_samples=8192;

post_Afout_Pfout(dx,dy,c,dt,x_ref,x_recorder,y_ref,y_recorder,x_bron,y_bron,recorder,recorder_ref,n_of_samples,tijdreeks,fc)
post_Afout_Pfout(dx,dy,c,dt,x_ref2,x_recorder2,y_ref2,y_recorder2,x_bron,y_bron,recorder2,recorder2_ref,n_of_samples,tijdreeks,fc)

rmpath(directory);