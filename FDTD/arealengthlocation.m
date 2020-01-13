function arealengthlocation(a,dx,dy,n)
a=0.19;
dx = 0.1;
dy = 0.1;
n = 6;
c = n/2+0.5;
global mx my area lxl lxr lyu lyd
mx=[];
my=[];
area=[];
lxl=[];
lxr=[];
lyu=[];
lyd=[];
for i = 1:n
    for j = 1:n
        x = (i-c)*dx;
        y = (c-j)*dy;
        x1=x-dx/2;
        x2=x+dx/2;
        y1=y-dy/2;
        y2=y+dy/2;
        if x1^2+y1^2 <= a^2 || x2^2+y1^2 <= a^2 || x1^2+y2^2 <= a^2 || x2^2+y2^2 <= a^2 
             mx = [mx i];
             my = [my j];
        end
    end
end
% for v=1:size(mx,2)
%     i = mx(1,v);
%     j = my(1,v);
%     x = (i-c)*dx;
%     y = (c-j)*dy;
%     x1=x-dx/2;
%     x2=x+dx/2;
%     y1=y-dy/2;
%     y2=y+dy/2;
%     s = dx/2000;
%     lx = [x1:s:x2];
%     ly = [y1:s:y2];
%     tt = dx/s+1;
%     opp = 0;
%     for ii=1:tt
%         for jj=1:tt
%             if lx(ii)^2+ly(jj)^2 <=a^2
%                 opp=opp+s^2;
%             end
%         end
%     end
%     sl = dx/1000000;
%     ttl = dx/sl+1;
%     l = [x1:sl:x2];
%     tll = 0;
%     tlr = 0;
%     tlu = 0;
%     tld = 0;
%     for iii=1:ttl
%             if x1^2+l(iii)^2 <=a^2
%                 tll = tll + sl;
%             end
%             if x2^2+l(iii)^2 <=a^2
%                 tlr = tlr + sl;
%             end
%             if l(iii)^2+y1^2 <=a^2
%                 tld = tld + sl;
%             end
%             if l(iii)^2+y2^2 <=a^2
%                 tlu = tlu + sl;
%             end
%     end 
%     area = [area a];
%     lxl = [lxl tll];
%     lxr = [lxr tlr];
%     lyu = [lyu tlu];
%     lyd = [lyd tld];
% end
mx
my
area
lxl


                    
            