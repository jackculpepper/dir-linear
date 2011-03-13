
loglike_ais = zeros(Btest,length(t_range));
loglike_ais_gauss = zeros(Btest,length(t_range));
loglike_act = zeros(Btest,length(t_range));

loglike_ais_mean = zeros(1,length(t_range));
loglike_ais_gauss_mean = zeros(1,length(t_range));
loglike_act_mean = zeros(1,length(t_range));

ratio = zeros(1,length(t_range));


%% load testing data
load_Xtest


t_start = tic();

for i = 1:length(t_range)
    %% number of intermediate distributions for AIS
    opts.T = ceil(t_range(i));

    for b = 1:Btest
        switch model
            case 'gauss'
                loglike_ais(b,i) = logZ(opts, Xtest(:,b)*ones(1,opts.BatchSize), phi, sigma);
                loglike_ais_gauss(b,i) = AIS_gauss(opts, Xtest(:,b)*ones(1,opts.BatchSize), phi, sigma);
            case 'laplace'
                loglike_ais(b,i) = logZ(opts, Xtest(:,b)*ones(1,opts.BatchSize), phi, sigma, lambda);
                loglike_ais_gauss(b,i) = AIS_gauss(opts, Xtest(:,b)*ones(1,opts.BatchSize), phi, sigma, lambda);
            case 'bins'
                opts.x0 = 2*rand( size(phi,2), opts.BatchSize )-1;
                loglike_ais(b,i) = logZ(opts, Xtest(:,b)*ones(1,opts.BatchSize), phi, sigma, gamma);
                loglike_ais_gauss(b,i) = AIS_gauss(opts, Xtest(:,b)*ones(1,opts.BatchSize), phi, sigma, gamma);
        end

        loglike_act(b,i) = logZ_gauss(Xtest(:,b), phi, sigma);

        t_c = toc(t_start);
        
        fprintf('T %07d Sample log likelihood via AIS: %f in %f sec\n', opts.T, loglike_ais(b, i), t_c);
        fprintf('T %07d Sample log likelihood (Gauss): %f\n', opts.T, loglike_act(b, i));
    end

    loglike_ais_mean(i) = mean(loglike_ais(:,i));
    loglike_ais_gauss_mean(i) = mean(loglike_ais_gauss(:,i));
    loglike_act_mean(i) = mean(loglike_act(:,i));

    fprintf('T %07d Average Log likelihood via AIS: %f in %f sec\n', opts.T, loglike_ais_mean(i), t_c);
    fprintf('T %07d Average Log likelihood via AIS gauss: %f in %f sec\n', opts.T, loglike_ais_gauss_mean(i), t_c);
    fprintf('T %07d Average Log likelihood (Gauss): %f\n', opts.T, loglike_act_mean(i));

    ratio(i) = loglike_ais_mean(i) / loglike_act_mean(i);
    fprintf('Ratio: %f\n', ratio(i));

    sfigure(14);
    semilogx(t_range(1:i), ratio(1:i), '.-');
    title('Ratio of estimated to true log likelihood vs. number of intermediate distributions');
    xlabel('Number of intermediate distributions');
    ylabel('Ratio of estimated to true log likelihood');

    sfigure(15);
    semilogx(t_range(1:i), loglike_ais_mean(1:i), '.-', t_range(1:i), loglike_ais_gauss_mean(1:i), '.-r', t_range(1:i), loglike_act_mean(1:i), '--');
    legend('HAIS', 'AIS Gauss', 'Gaussian', 'Location', 'Best');
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


