function g = dEdX_laplace(x0, X, phi, sigma, lambda);

[L M] = size(phi);

B = size(X,2);
a = reshape(x0,M,B);

E = X - phi*a;

g = -phi'*inv(sigma)*E + lambda*sign(a);

