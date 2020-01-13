function piberekenen
%     lx = -1 + (2)*rand(100000000,1);
%     ly = -1 + (2)*rand(100000000,1);
%     teller=0;
%     for i=1:100000000
%     if lx(i)^2+ly(i)^2 <= 1^2
%         teller=teller+1;
%     end
%     end
%     (teller)/100000000*2^2/1^2
%     end
s = 0.0005;
lx = [-1:s:1];
ly = [-1:s:1];
tt = 2/s+1;
a = 0;
for i=1:tt
    for j=1:tt
        if lx(i)^2+ly(j)^2 <=1^2
            a=a+s^2;
        end
    end
end
vpa(a)