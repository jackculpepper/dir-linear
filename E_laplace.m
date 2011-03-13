function f = E_laplace(x0, X, phi, sigma, lambda);

[L M] = size(phi);
B = size(X,2);

a = reshape(x0,M,B);

%% likelihood normalization factor
z_l = log( (2*pi)^(L/2) * det(sigma)^(1/2) );
%% prior normalization factor
z_p = M*log(2/lambda);

E = X - phi*a;

%% likelihood energy
f_l = 0.5 * sum( E .*(inv(sigma) * E), 1 );
%% prior energy
f_p = lambda * sum(abs(a));

f = f_l + z_l + f_p + z_p;
f = f(:)';

