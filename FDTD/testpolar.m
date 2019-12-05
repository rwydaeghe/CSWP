c=340;
dr=0.1;
Nr=10;
roff=0.0001;
ro=roff+Nr*dr;
dh=0.1;
Nh=round(2*pi/dh);
dh=2*pi/Nh;

global or ot p
p = zeros(Nr, Nh); %low Nr index (row) low r low Nh index (row) low h
or = zeros(Nr+1, Nh);
ot = zeros(Nr,Nh);

CFL=1; %Courant getal - Courant number

dt=CFL/c*(1/dr^2+(2/dr/dh)^2)^(-1/2); %tijdstap - time step

nt=500/CFL; %aantal tijdstappen in simulatie - number of time steps

A=1;fc=200;t0=2.5E-10;sigma=5E-10;

%film
%movie
mov=moviein(nt);

bront = zeros(nt,1);
tijdreeks=zeros(nt,1);
bron=0;

figure;

R=[roff+dr/2:dr:ro-dr/2];
T=[dh/2:dh:2*pi-dh/2]*180/pi;

%TIJDSITERATIE------------------------------------------------------
%TIME ITTERATION----------------------------------------------------
for it=1:nt
   t = (it-1)*dt; tijdreeks(it)=t;
	disp([num2str(it) '/' num2str(nt)]);
  
    bron=A*cos(2*pi*fc*(t-t0))*exp(-((t-t0)^2)/(sigma)) %bron updaten bij nieuw tijd - update source for new time
   for i=1:Nh
       or(1,i)=bron;
   end
   %p(x_bron,y_bron) = p(x_bron,y_bron)+bron; %druk toevoegen bij de drukvergelijking op bronlocatie - adding source term to propagation
    %ot = zeros(Nr,Nh);
   newstep(Nr,Nh,c,dr,dh,dt,bron,roff);   %propagatie over 1 tijdstap - propagate over one time step
    %ot = zeros(Nr,Nh);

%   p(isnan(p))=0;
%    or(isnan(or))=0;
%    ot(isnan(ot))=0;
   %voorstellen drukveld
   %presenting the p field
   polarPcolor(R,T,p); shading interp; caxis([-0.02*A 0.02*A]); title([num2str(it) '/' num2str(nt)]);
  
   %mov(it) = getframe; %wegcommentarieren voor simulatie vlugger te laten lopen - if this line is removed simulation runs much faster
end
movie(mov)







