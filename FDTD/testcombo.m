clear all

c=340; %speed of sound
dx=0.1; %spatial discretisation step size x direction
dy=dx; %spatial discretisation step size y direction (should be same as x)
a=1; %radius of the circle
A=1; %amplitude plane wave
k=20; %k plane wave

nx=round(6*a/dx);
if mod(nx,2)==1
    nx = nx+1; %number of cells in x direction
end

ny=nx; %number of cells in y direction (same as nx for practical reasons)

nc=round(a/dx)+1;
if mod(nc,2)==1
    nc = nc+1; %number of cells to make square around circle
end

CFL=1; %Courant number
dt=CFL/(c*sqrt((1/dx^2)+(1/dy^2))); %time step
nt=400/CFL; %number of time steps

%possibility to be done: interpolate for exact point

x_c=round(nx/2); y_c=round(ny/2);



x_recorder=x_c; y_recorder=y_c+round(2*a/dx); %location receiver 1
x_ref=x_c;y_ref=y_c+round(2.5*a/dx); %location reference receiver 1

x_recorder2=x_c; y_recorder2=y_c-round(2*a/dx); %location receiver 2
x_ref2=x_c;y_ref2=y_c-round(2.5*a/dx); %location reference receiver 2

x_recorder3=x_c-round(2*a/dx); y_recorder3=y_c; %location receiver 3
x_ref3=x_c-round(2.5*a/dx);y_ref3=y_c; %location reference receiver 3

x_recorder4=x_c+round(2*a/dx); y_recorder4=y_c; %location receiver 4
x_ref4=x_c+round(2.5*a/dx);y_ref4=y_c; %location reference receiver 4


