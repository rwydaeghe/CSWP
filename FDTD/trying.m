global mx my
q = ones(26,26)
for v=1:size(mx,2)
    q(my(v),mx(v)) = 0;
end
gridplot(q)
    