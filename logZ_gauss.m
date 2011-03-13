function f = logZ_gauss(X, phi, sigma, varargin);

[L M] = size(phi);
sigma_posterior = phi*phi' + sigma;
z = log( (2*pi)^(L/2) * det(sigma_posterior)^(1/2) );
f = trace(-0.5*X'*inv(sigma_posterior)*X - z);

