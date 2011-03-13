
phi = randn(L,M);
phi = 10*phi*diag(1./sqrt(sum(phi.^2)));

load_X


%% init the latent vars to the prior
switch model
    case 'laplace'
        a = rand(M, B);
        a = sign(randn(M,B)) .* expinv(a, 1/lambda);
    case 'gauss'
        a = randn(M, B);
end


p = zeros(M, B);

update = 1;

p_log = [];
eta_log = [];
stepadj_log = [];
energies_log = [];
updatelength_log = [];
loglike_ais_log = [];


