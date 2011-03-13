
M = 5;
L = 3;
B = 2;

X = randn(L,B);
phi = randn(L,M);
a = randn(M,B);

sigma = randn(L)/sqrt(L);
sigma = sigma'*sigma;

tic
checkgrad('E_dEdTheta', phi(:), 1e-4, X, a, sigma)
toc

