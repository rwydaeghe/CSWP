function s = analytic(E0,a,k,n,r,phi)
    c=340;
    s=0;
    for l = 0:n
        v=2;
        if l == 0
            v=1;
        end
         a1 = (besselj(l, k*a)*l)/k/a - besselj(1 + l, k*a);
         a2 = (besselh(l, 1, k*a)*l)/k/a - besselh(1 + l, 1, k*a);
        s = s + E0*v*(1i)^(l)*(1*besselj(l,k*r)-1*a1/a2*besselh(l,1,k*r))*cos(l*(phi-pi));
    end
end