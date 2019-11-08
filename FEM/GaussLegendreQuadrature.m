% Syntax: [W,X] = GaussLegendreQuadrature(n)
% This function returns an n-point Gauss-Legendre quadrature rule that exactly 
% integrates polynomials up to degree 2n-1.
% Output: X = the quadrature nodes
%         W = the quadrature weights
% These rules are for integration over the interval x in [-1,1].
% Rules are available for any positive and integer n.
% Example: 
% [W,X] = GaussLegendreQuadrature(6)
% Int = sum(W.*X.^10)
% which should be close to the analytical result 2/11.
% int(x^10,x=-1..1) = 2/11
function [W,X] = GaussLegendreQuadrature(n)
NofAppNodes = floor(n/2);
X = EstimateLegendreNodes(n);

% Refine approximated nodes with Newton iteration
Tol = 1e-14;
for k = 1:NofAppNodes
    [x, Update] = UpdateNode(n, X(k));
    while(abs(Update) > Tol)
        [x, Update] = UpdateNode(n, x);
    end
    [x, Update] = UpdateNode(n, x);
    X(k) = x;
end

% Now calculate the weights of the strictly positive nodes
P = LegendreP_l(n+1, X);
W = 2*(1-X.^2)/(n+1)^2./P.^2;

% And use symmetry to expand
if(mod(n,2) == 1) % if odd
    X0 = 0;
    W0 = 2/(n+1)^2./LegendreP_l(n+1, X0).^2;
else % if even
    X0 = [];
    W0 = [];
end
W = [+W,W0,flipdim(W,2)].';
X = [-X,X0,flipdim(X,2)].';
return

function Nodes = EstimateLegendreNodes(N)
% This function gives good estimates for the Legendre zeros
NofNodes =  floor(N/2);
Nodes = (N*N+2*N-0.478683) / (N*N+2*N-0.341487) * cos((4*[1:NofNodes]-1)*pi/(4*N+2));
return

function [x, Update] = UpdateNode(n, x)
Leg = LegendreP_l(n,x);
Update = Leg*(1-x*x)/(n+1) / (x*Leg - LegendreP_l(n+1,x));
x = x - Update;
return

function P = LegendreP_l(l, x)
% This function calculates the Legendre polynomial of degree l for the 
% arguments contained in the vector x.
P = ones(size(x));
Poud = zeros(size(x));
for n = 0:(l-1)
    Pzeeroud = Poud;
    Poud = P;
    P = (1/(n+1)) * (-n*Pzeeroud + (2*n+1) * x .* Poud);
end
return