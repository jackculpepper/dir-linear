function g = dEdTheta_linear(x0, X, a, sigma)

M = size(a,1);
[L B] = size(X);

phi = reshape(x0, L, M);

E = X - phi*a;

g = -inv(sigma)*E*a';

g = g(:);


