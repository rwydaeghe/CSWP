function reference(dx,dt,n,nt,A,k,c)
global oxref pref
for it=1:nt+200
   t = (it-1)*dt; tijdreeks(it)=t;
	disp([num2str(it) '/' num2str(nt)]);
  
    bron=A*sin(k*c*(t)); %bron updaten bij nieuw tijd - update source for new time
   
    oxref(it,1) = oxref(it,1)+bron; %druk toevoegen bij de drukvergelijking op bronlocatie - adding source term to propagation
    for i = 1:n
        pref(it+1,i) = pref(it,i) - c^2*dt/dx*(oxref(it,i+1)-oxref(it,i));
    end
    for i = 2:n
        oxref(it+1,i) = oxref(it,i) - dt/dx*(pref(it+1,i)-pref(it+1,i-1));
    end
    oxref(it+1,n+1) = 1/(1+c*dt/dx)*((1-c*dt/dx)*oxref(it,n+1)+2*dt/dx*pref(it+1,n));
%     plot(1:dx:1+n*dx,oxref(it+1,:))
%    mov(it) = getframe;
end
