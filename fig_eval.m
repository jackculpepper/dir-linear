
loglike_ais = zeros(Btest,length(t_range));
loglike_gauss = zeros(Btest,length(t_range));

loglike_ais_mean = zeros(1,length(t_range));
loglike_gauss_mean = zeros(1,length(t_range));

ratio = zeros(1,length(t_range));


%% load testing data
load_Xtest


t_s = tic();

for i = 1:length(t_range)
    %% number of intermediate distributions for AIS
    opts_hais.T = ceil(t_range(i));

    for b = 1:Btest
        %% copy test point b once for each particle
        Xin = Xtest(:,b)*ones(1,opts_hais.B);

        t_b = tic();
        switch model
            case 'gauss'
                loglike_ais(b,i) = hais(opts_hais, Xin, phi, sigma);
            case 'laplace'
                loglike_ais(b,i) = hais(opts_hais, Xin, phi, sigma, lambda);
        end

        loglike_gauss(b,i) = logZ_gauss(Xtest(:,b), phi, sigma);

        fprintf('T %07d %d/%d Sample log likelihood: %f (%f) in %.2fs (%.2fs total)\n', ...
            opts_hais.T, b, Btest, loglike_ais(b, i), loglike_gauss(b, i), toc(t_b), toc(t_s));
    end

    loglike_ais_mean(i) = mean(loglike_ais(:,i));
    loglike_gauss_mean(i) = mean(loglike_gauss(:,i));

    fprintf('T %07d Average Log likelihood via AIS: %f (%f) in %f sec\n', ...
        opts_hais.T, loglike_ais_mean(i), loglike_gauss_mean(i), t_c);

    ratio(i) = loglike_ais_mean(i) / loglike_gauss_mean(i);
    fprintf('Ratio: %f\n', ratio(i));

    sfigure(14);
    semilogx(t_range(1:i), ratio(1:i), '.-');
    title('Ratio of estimated to true log likelihood vs. number of intermediate distributions');
    xlabel('Number of intermediate distributions');
    ylabel('Ratio of estimated to true log likelihood');

    sfigure(15);
    semilogx(t_range(1:i), loglike_ais_mean(1:i), '.-', t_range(1:i), loglike_gauss_mean(1:i), '--');
    legend('HAIS', 'Gaussian', 'Location', 'Best');
    title('Estimated average log likelihood vs. number of intermediate distributions');
    xlabel('Number of intermediate distributions');
    ylabel('Log likelihood');
    axis tight;

    sfigure(16);
    semilogx(t_range(1:i), loglike_ais(:,1:i), '.-');
    title('Estimated log likelihood vs. number of intermediate distributions');
    xlabel('Number of intermediate distributions');
    ylabel('Log likelihood');
    axis tight;

    drawnow;
end


eval(sprintf('save state/%s/matlab_up=%06d_ais.mat', paramstr, update)); 


