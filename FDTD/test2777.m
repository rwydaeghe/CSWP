global q
q = ones(160,160)
circlewithpolar(0.1,0.1,1,160);
c = n/2+0.5;
for j=1:160-1
    for i=1:160
        if q(j,i) ~= q(j+1,i)
            q(j,i) = (c-j)*dy-dy/2;
%             q(j,i) = (c-j)*dy
        end
    end
end
% circlewithpolar(dx,dy,a,n);