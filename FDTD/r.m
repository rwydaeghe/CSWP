c=340;
dr=0.1;
ro=0.01+20*dr;
Nr=20;
dh=0.1/ro;
Nh=round(2*pi/dh);
dh=2*pi/Nh;

global or ot p
p = zeros(Nr, Nh); %low Nr index (row) low r low Nh index (row) low h
or = zeros(Nr+1, Nh);
ot = zeros(Nr,Nh);

%film
%movie
mov=moviein(nt);

CFL=1; %Courant getal - Courant number

dt=CFL/(c*sqrt((1/dr^2)+(1/(0.01^2*dh^2))); %tijdstap - time step

nt=400/CFL; %aantal tijdstappen in simulatie - number of time steps







