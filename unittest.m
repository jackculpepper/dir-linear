
clear

%% this many pca components
L = 25;
L = 36;
L = 128;


M = 144;

%% this many pixels
J = 64;
J = 256;
Jsz = sqrt(J);

Mrows = sqrt(L);
Mrows = sqrt(M);


display_every = 50;
display_every = 100;
display_every = 10;
%display_every = 1;

save_every = display_every;
test_every = 200000;
test_every = 100000;

reject_test = 0;

buff = 4;



noise_var = 0.1;
noise_var = 0.01;
%noise_var = 0.001;

sigma = sqrt(noise_var) * eye(L);
sigma = noise_var*eye(L);

lambda = 1;


Btest = 100;


datasource = 'movies';
datasource = 'images';
datasource = 'vanHat';

%% how many image patches to take from each vanHateren image. it takes time
%% to load a new image file; increasing this vastly speeds up the loading proc
patches_per_image_test = 1;
patches_per_image = 100;

learn_sigma = 0;
log_images = 1;


model = 'gauss';
model = 'laplace';


paramstr = sprintf('L=%03d_M=%03d_%s',L,M,datestr(now,30))

%% lbfgs options
opts_lbfgs = lbfgs_options('iprint', -1, 'maxits', 10, ...
                           'factr', 1e-1, 'cb', @cb);


%% hais options
opts_hais.BatchSize = 200;		% number of particles
opts_hais.DataSize = M;

switch model
    case 'gauss'
        opts_hais.E = @E_gauss;
        opts_hais.dEdX = @dEdX_gauss; 
    case 'laplace'
        opts_hais.E = @E_laplace;
        opts_hais.dEdX = @dEdX_laplace; 
end

opts_hais.epsilon = 0.1;
opts_hais.epsilon = 0.01;
opts_hais.Debug = 0;
opts_hais.sumlogZ = 1;
opts_hais.interpType = 'cos';
opts_hais.interpType = 'lin';

 
%% langevin options
langevin_epsilon = 0.01;         % take simulation steps of length epsilon
langevin_nsteps = 25;            % steps to take between every update of phi
langevin_nsteps = 10;            % steps to take between every update of phi
langevin_alpha = 0.5;
langevin_beta = 1 - exp(langevin_epsilon * log(langevin_alpha) );

B = 10000;
B = 50000;

reinit

t_range = 10.^5;

num_trials = 1000;

flag_reject = 0;
for i = 1:10 ; learn ; end

if reject_test
    flag_reject = 1;
    for i = 1:2 ; learn ; end
end

t_range = 10.^[2:0.5:5];

fig_ais

