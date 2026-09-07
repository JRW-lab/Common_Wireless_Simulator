function I = kink_split_integral(fun, a, b, kinks, n)
% Composite n-point Gauss-Legendre quadrature of fun over [a,b], splitting
% the domain at every point in `kinks` that lies strictly inside (a,b), so
% each sub-interval handed to Gauss-Legendre is smooth. This gives
% spectral (near machine-precision) accuracy per sub-interval instead of
% the slow, first-order convergence a naive Riemann sum or a single global
% Gauss-Legendre application would have across a kink (a discontinuity in
% value or derivative, e.g. a pulse's truncation edge).
%
% fun must accept a column vector of evaluation points and return a column
% vector of the same size (vectorized).

if a > b
    tmp = a; a = b; b = tmp;
end

pts = unique([a; b; kinks(:)]);
pts = pts(pts >= a - 1e-12 & pts <= b + 1e-12);
pts = sort(pts);
pts([false; diff(pts) < 1e-14]) = [];

if numel(pts) < 2
    I = 0;
    return;
end

[nodes,weights] = gauss_legendre_nodes(n);

I = 0;
for i = 1:numel(pts)-1
    lo = pts(i);
    hi = pts(i+1);
    if hi <= lo
        continue;
    end
    xq = (hi-lo)/2 * nodes + (hi+lo)/2;
    fq = fun(xq);
    I = I + (hi-lo)/2 * (weights.' * fq(:));
end

end
