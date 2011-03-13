function [f, g] = wrap_dEdTheta(x0, X, a, sigma);

[L, B] = size(X);
M = size(a, 1);

phi = reshape(x0, L, M);

f = E_gauss(a(:), X, phi, sigma);
g = dEdTheta(x0, X, a, sigma);

f = sum(f(:));
g = g(:);

