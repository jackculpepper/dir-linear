function g = dEdX_gauss(x0, X, phi, sigma);

[L B] = size(X);
M = size(phi,2);

a = reshape(x0, M, B);

E = X - phi*a;

g = -phi'*inv(sigma)*E + a;