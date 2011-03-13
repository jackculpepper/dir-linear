%% inverse cdf for the exponential distribution
function x = expinv(p,mu)

k = (0 < p & p < 1);

if all(k(:))
    q = -log(1-p);
else
    q = zeros(size(p));
    q(k) = -log(1-p(k));
    q(p == 1) = Inf;
    q(p < 0 | 1 < p | isnan(p)) = NaN;
end

mu(mu <= 0) = NaN;
x = mu .* q;


