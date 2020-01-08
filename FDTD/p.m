% global q
% c = n/2+0.5;
% for j=1:160
%     for i=1:160-1
%         if q(j,i) ~= q(j,i+1)
%             x = (i-c)*dx+dx/2
%             y = (c-j)*dy
%         end
%     end
% end
circlewithpolar(dx,dy,a,n);