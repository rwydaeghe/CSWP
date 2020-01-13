function post_Afout_Pfout(n,a,dx,dy,c,dt,x_ref,x_recorder,y_ref,y_recorder,x_bron,y_bron,recorder,recorder_ref,n_of_samples,tijdreeks,k)
fc=k*c/2/pi;
%generatie frequentie-as - generate frequency axes
maxf=20*c/2/pi;
minf=c/2/pi;
df=(maxf-minf)/n_of_samples;
fas=zeros(n_of_samples,1);
for loper=1:n_of_samples
   fas(loper,1)=df*(loper-1)+minf;
end

%amplitudeverhouding en faseverhouding analytisch - analytical amplitude
%ratio and phase difference
cc = n/2+0.5;
xref = (x_ref-cc)*dx;
yref = (cc-y_ref)*dy;
xrecorder = (x_recorder-cc)*dx;
yrecorder = (cc-y_recorder)*dy;
r1=sqrt((xref)^2+(yref)^2);
theta1=atan2(yref,xref);
theta1 = theta1.*(theta1 >= 0) + (theta1 + 2*pi).*(theta1 < 0);
r2=sqrt((xrecorder)^2+(yrecorder)^2);
theta2=atan2(yrecorder,xrecorder);
theta2 = theta2.*(theta2 >= 0) + (theta2 + 2*pi).*(theta2 < 0);
aantalcellengepropageerd=sqrt((x_recorder-x_ref)^2+(y_recorder-y_ref)^2);

k=2*pi*fas/c;
f1=zeros(n_of_samples,1);
f2=zeros(n_of_samples,1);
for i=1:n_of_samples
    f1(i)=analytic(1,a,k(i),50,r1,theta1);
    f2(i)=analytic(1,a,k(i),50,r2,theta2);
end
Averhouding_theorie=abs(f1./f2);
Pverschil_theorie=(angle(f1))-(angle(f2));

%amplitudeverhouding en faseverhouding FDTD - amplitude ratio and phase
%difference from FDTD
fftrecorder=fft(recorder,n_of_samples);
fftrecorder_ref=fft(recorder_ref,n_of_samples);
Averhouding_FDTD=abs(fftrecorder_ref./fftrecorder);
Pverschil_FDTD=(angle(fftrecorder_ref))-(angle(fftrecorder));

% figure;subplot(2,3,1);plot(tijdreeks,recorder);title('t recorder');subplot(2,3,2);plot(fas,abs(fftrecorder));title('fft recorder abs');xlim([0.5*fc 1.5*fc]);
% subplot(2,3,3);plot(fas,unwrap(angle(fftrecorder)));title('fft recorder phase');xlim([0.5*fc 1.5*fc]);
% subplot(2,3,4);plot(tijdreeks,recorder_ref);title('t recorder ref');subplot(2,3,5);plot(fas,abs(fftrecorder_ref));title('fft recorder ref abs');xlim([0.5*fc 1.5*fc]);
% subplot(2,3,6);plot(fas,unwrap(angle(fftrecorder_ref)));title('fft recorder phase');xlim([0.5*fc 1.5*fc]);

%vergelijking analytisch-FDTD - comparison analytical versus FDTD
ka=a*2*pi*fas/c;



 Averhoudingrel=abs(Averhouding_FDTD./Averhouding_theorie);
 Averhouding=1+((Averhoudingrel-1)/aantalcellengepropageerd);

figure;subplot(2,1,1);plot(ka,Averhouding);xlim([1 20]);title('Amplitudeverhouding FDTD/analyt. per ka');ylabel('verhouding');xlabel('ka');ylim([0.9 1.5]);

%Pverhouding=Pverschil_FDTD-Pverschil_theorie;
Pverhouding=(Pverschil_FDTD+Pverschil_theorie)/aantalcellengepropageerd;
subplot(2,1,2);plot(ka,Pverhouding);xlim([1 20]);title('Faseverschil FDTD/analyt. per ka');ylabel('verschil');xlabel('ka');ylim([-0.1 0.1]);
% ka(1000:1100)
% abs(f2(1000:1100))
% abs(fftrecorder(1000:1100))
%  angle(f2(1000:1100)./f1(1000:1100))
%  angle(fftrecorder(1000:1100)./fftrecorder_ref(1000:1100))