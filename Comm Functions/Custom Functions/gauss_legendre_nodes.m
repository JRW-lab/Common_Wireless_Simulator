function [nodes,weights] = gauss_legendre_nodes(n)
% Nodes and weights of n-point Gauss-Legendre quadrature on [-1,1], via the
% Golub-Welsch algorithm (eigendecomposition of the Jacobi matrix built
% from the three-term recurrence of the Legendre polynomials).
%
% Standard, well-established numerical method (e.g. Golub & Welsch 1969).

beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^-2);
J = diag(beta,1) + diag(beta,-1);
[V,D] = eig(J);
[nodes,idx] = sort(diag(D));
weights = 2 * (V(1,idx).^2);
weights = weights(:);

end
