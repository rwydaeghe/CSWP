% Syntax: [W,X] = GaussLegendreQuadrature01(n)
% This function returns an n-point Gauss-Legendre quadrature rule that exactly 
% integrates polynomials up to degree 2n-1.
% Output: X = the quadrature nodes
%         W = the quadrature weights
% These rules are for integration over the interval x in [0,1].
% Rules are available for any positive and integer n.
% Example: 
% [W,X] = GaussLegendreQuadrature(6)
% Int = sum(W.*X.^10)
% which should be close to the analytical result 1/11.
% int(x^10,x=0..1) = 1/11
function [W,X] = GaussLegendreQuadrature01(n)
[W,X] = GaussLegendreQuadrature(n);
W = W/2;
X = (X+1)/2;
return