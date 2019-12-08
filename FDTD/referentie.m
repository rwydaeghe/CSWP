nx = 20;
nsf = 20;
n = nx;
o = zeros(1,n+1);
p = zeros(1,n)
dx = 0.4;
c = 340;
CFL=1; %Courant getal - Courant number

dt=CFL/(c*sqrt((1/dx^2)+(1/dx^2))); %tijdstap - time step

nt=400;
A=1; k=1;


for it=1:nt-1
   t = (it-1)*dt; tijdreeks(it)=t;
	disp([num2str(it) '/' num2str(nt)]);
  
    p(1,1:n) = p(1,1:n) - c^2*dt/dx*(o(1,2:n+1)-o(1,1:n));
    o(1,1) = A*sin(k*c*t)
    o(1,2:n) = o(1,2:n) - dt/dx*(p(1,2:n)-p(1,1:n-1));
    o(1,n+1) = 1/(1+c*dt/dx)*((1-c*dt/dx)*o(1,n+1)+2*dt/dx*p(1,n));
    plot(1:dx:1+n*dx,o)
   mov(it) = getframe;
end
    
    
    
    
    