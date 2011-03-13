function logZ = hais(opts, theta, varargin);

ET = opts.E;
dETdX = opts.dEdX;


%% fix the proposal distribution to be univariate gaussian
E0 = @E_gauss_unit_norm;
dE0dX = @dEdX_gauss_unit_norm;
logZ0 = @logZ_gauss_unit_norm;


B = opts.B;
T = opts.T;
M = opts.M;

%% fix the parameters for the sampler
epsilon = 0.01;
alpha = 0.5;
beta = 1 - exp(log(alpha) * epsilon);


%% initialize position and momentum
x0 = randn(M, B);
p0 = randn(M, B);

x = x0;
p = p0;
logw = 0;

num_rej = 0;
num_tot = 0;


%% evaluate the target and proposal energy functions
E0t = E0(x, theta, varargin{:});
ETt = ET(x, theta, varargin{:});

Em0 = E0t;

for t = 2:T

    %% mix linearly between the target and proposal
    gamma = t/T;
    Em1 = (1 - gamma)*E0t + gamma*ETt;

    logw = logw - Em1 + Em0;

    %% partial momentum refresh
    p = -sqrt(1 - beta)*p + sqrt(beta)*randn(M, B);



    %% langevin dynamics -- transition

    %% leapfrog 1/2 step
    p0 = p;
    x0 = x;
    E0t0 = E0t;
    ETt0 = ETt;
    x = x + epsilon/2 * p;


    dEm0dX = dE0dX(x, theta, varargin{:});
    dEmTdX = dETdX(x, theta, varargin{:});

    %% mix target and proposal gradients
    dE = (1 - gamma)*dEm0dX + gamma*dEmTdX;

    %% leapfrog 1/2 step
    p = p - epsilon * dE;  
    x = x + epsilon/2 * p;


    %% reflect momentum
    p = -p;


    E0t = E0(x, theta, varargin{:});
    ETt = ET(x, theta, varargin{:});
    Em0 = (1 - gamma)*E0t + gamma*ETt;


    %% MH rejection
    E_delta = E_gauss_unit_norm(p) + Em0 - E_gauss_unit_norm(p0) - Em1;

    p_acc = exp(-E_delta);
    idx = p_acc < rand(1, B);

    p(:,idx) = p0(:,idx);
    x(:,idx) = x0(:,idx);

    Em0(idx) = Em1(idx);
    E0t(idx) = E0t0(idx);
    ETt(idx) = ETt0(idx);

    num_rej = num_rej + sum(idx);
    num_tot = num_tot + size(p,2);

    if opts.flag_debug ; fprintf('\rt %d / %d', t, T); end
end

logZweights = logZ0(M, varargin{:}) + logw;
k = max(logZweights);
logZ = log(sum(exp(logZweights - k))/B) + k;

if isinf(logZ)
    keyboard
end

if opts.flag_debug
    fprintf(' reject fraction %f\n', num_rej / num_tot);
end




%% univariate normal logZ, E, and dEdX

function dE = dEdX_gauss_unit_norm(X, varargin)
dE = X;

function E = E_gauss_unit_norm(X, varargin)
E = 0.5*sum(X.^2, 1);

function logZ = logZ_gauss_unit_norm(M, varargin )
logZ = (M/2) * log( 2*pi );


