function [f, g] = E_dEdX_gauss(x0, X, phi, sigma);

g = dEdX_gauss(x0, X, phi, sigma);
f = E_gauss(x0, X, phi, sigma);

f = sum(f(:));
g = g(:);

