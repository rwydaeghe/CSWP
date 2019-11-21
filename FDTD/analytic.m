function s = analytic(E0,a,k,n,r,phi)
    c=340;
    s=0;
    for i = 0:n
        v=2;
        if i == 0
            v=1;
        end
        s = s + E0*v*(-1i)^(i)*(besselj(i,k*r)-besselj(i,k*a)/besselh(i,2,k*a)*besselh(i,2,k*r))*cos(n*phi);
    end
end