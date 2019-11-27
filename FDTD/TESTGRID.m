R=6000; %radius
S=15;   %num circ.lines
N=16;   %num ang.lines
sect_width=2*pi/N;    
offset_angle=(deg2rad(11.25)):sect_width:2*pi-sect_width+deg2rad(11.25);
%------------------
r=linspace(4000,R,S+1);
w=0:.01:2*pi;
clf %remove if needed
figure (1);
hold on
axis equal
for n=2:length(r)
      plot(real(r(n)*exp(1i*w)),imag(r(n)*exp(1i*w)),'k--')
end 
% 
for n=1:length(offset_angle)
      plot(real([4000 R]*exp(1i*offset_angle(n))),imag([4000 R]*exp(1i*offset_angle(n))),'k-')
end