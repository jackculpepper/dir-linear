function f = E_gauss(x0, X, phi, sigma);

[L B] = size(X);
M = size(phi,2);

a = reshape(x0, M, B);

%% likelihood normalization factor
z_l = log( (2*pi)^(L/2) * det(sigma)^(1/2) );
%% prior normalization factor
z_p = log( (2*pi)^(M/2) );

E = X - phi*a;

%% likelihood energy
f_l = 0.5 * sum( E .*(inv(sigma) * E), 1 );
%% prior energy
f_p = 0.5 * sum( a.^2, 1 );

f = f_l + z_l + f_p + z_p;
f = f(:)';

