
M = 5;
L = 3;
B = 2;

X = randn(L,B);
phi = randn(L,M);
a = randn(M,B);

sigma = 0.01*eye(L);

tic
checkgrad('E_dEdX_gauss', a(:), 1e-4, X, phi, sigma)
toc

