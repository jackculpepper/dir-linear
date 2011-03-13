function [f, g] = wrap_dEdX_laplace(x0, X, phi, sigma, lambda);

g = dEdX_laplace(x0, X, phi, sigma, lambda);
f = E_laplace(x0, X, phi, sigma, lambda);

f = sum(f(:));
g = g(:);

