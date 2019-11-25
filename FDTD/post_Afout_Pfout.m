function post_Afout_Pfout(a,dx,dy,c,dt,x_ref,x_recorder,y_ref,y_recorder,x_bron,y_bron,recorder,recorder_ref,n_of_samples,tijdreeks,k)
fc=k*c/2/pi;
%generatie frequentie-as - generate frequency axes
maxf=1/dt;
df=maxf/n_of_samples;
fas=zeros(n_of_samples,1);
for loper=1:n_of_samples
   fas(loper,1)=df*(loper);
end

%amplitudeverhouding en faseverhouding analytisch - analytical amplitude
%ratio and phase difference
r1=sqrt(((x_ref-x_bron)*dx)^2+((y_ref-y_bron)*dy)^2)
theta1=atan(((y_ref-y_bron)*dy)/((x_ref-x_bron)*dx))
r2=sqrt(((x_recorder-x_bron)*dx)^2+((y_recorder-y_bron)*dy)^2)
theta2=atan(((y_recorder-y_bron)*dy)/((x_recorder-x_bron)*dx))
aantalcellengepropageerd=sqrt((x_recorder-x_ref)^2+(y_recorder-y_ref)^2)

k=2*pi*fas/c;
f1 = analytic(1,a,k,10,r1,theta1)
f2 = analytic(1,a,k,10,r2,theta2)
Averhouding_theorie=abs(f1./f2);
Pverschil_theorie=unwrap(angle(f1))-unwrap(angle(f2));

%amplitudeverhouding en faseverhouding FDTD - amplitude ratio and phase
%difference from FDTD
fftrecorder=fft(recorder,n_of_samples);
fftrecorder_ref=fft(recorder_ref,n_of_samples);
Averhouding_FDTD=abs(fftrecorder_ref./fftrecorder);
Pverschil_FDTD=unwrap(angle(fftrecorder_ref))-unwrap(angle(fftrecorder));

figure;subplot(2,3,1);plot(tijdreeks,recorder);title('t recorder');subplot(2,3,2);plot(fas,abs(fftrecorder));title('fft recorder abs');xlim([0.5*fc 1.5*fc]);
subplot(2,3,3);plot(fas,unwrap(angle(fftrecorder)));title('fft recorder phase');xlim([0.5*fc 1.5*fc]);
subplot(2,3,4);plot(tijdreeks,recorder_ref);title('t recorder ref');subplot(2,3,5);plot(fas,abs(fftrecorder_ref));title('fft recorder ref abs');xlim([0.5*fc 1.5*fc]);
subplot(2,3,6);plot(fas,unwrap(angle(fftrecorder_ref)));title('fft recorder phase');xlim([0.5*fc 1.5*fc]);

%vergelijking analytisch-FDTD - comparison analytical versus FDTD
lambdaoverdx=(c./fas)./dx;
Averhoudingrel=(Averhouding_FDTD./Averhouding_theorie);
Averhouding=1+((Averhoudingrel-1)/aantalcellengepropageerd);

figure;subplot(2,1,1);plot(lambdaoverdx,Averhouding);xlim([5 20]);title('Amplitudeverhouding FDTD/analyt. per cel');ylabel('verhouding');xlabel('aantal cellen per golflengte');ylim([0.99 1.01]);
Pverhouding=(Pverschil_FDTD+Pverschil_theorie)/aantalcellengepropageerd;
subplot(2,1,2);plot(lambdaoverdx,Pverhouding);xlim([5 20]);title('Faseverschil FDTD/analyt. per cel');ylabel('verschil');xlabel('aantal cellen per golflengte');ylim([-0.01 0.01]);
