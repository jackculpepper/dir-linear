
for t = 1:num_trials

    if update == 1
        %% do many integration steps at the beginning, to mix the markov chain
        lst = 10 * langevin_nsteps;
    else
        %% randomize the number of integration steps
        lst = langevin_nsteps + sign(randn)*ceil(langevin_nsteps/2*rand);
    end

    tic

    if ~flag_reject
        %% rejectionless dynamics
        for i = 1:lst

            %% evaluate the gradient of the energy function wrt position vars
            switch model
                case 'laplace'
                    dEda = dEdX_laplace(a(:), X, phi, sigma, lambda);
                case 'gauss'
                    dEda = dEdX_gauss(a(:), X, phi, sigma);
            end

            %% draw from univariate gaussian to partially corrupt momentum
            R = randn(M, B);

            %% corrupt momentum
            p = sqrt(1-langevin_beta)*p + sqrt(langevin_beta)*R - ...
                langevin_epsilon*dEda;

            %% do not let the momentum grow too large
            idx = find( abs(p(:)) > 10 );
            p(idx) = sign(p(idx))*10;

            %% run dynamics forward in time
            a = a + langevin_epsilon*p;

            fprintf('\rlangevin step %d / %d', i, lst);

            if 1
                %% spend a few cpu cycles monitoring things during sampling
                E = X - phi*a;
                snr = 10 * log10(sum(X(:).^2) / sum(E(:).^2));
                fprintf(' snr %.4f', snr);
                fprintf(' a_abs %.4f', mean(abs(a(:))));
                fprintf(' p_var %.4f', var(p(:)));
            end
        end
    else
        %% regular dynamics with Metropolis-Hastings rejection criteria
        for i = 1:lst
            %% hamiltonian dynamics
            ainit = a;
            pinit = p;

            %% initial half step
            a = a + (1/2)*langevin_epsilon*p;

            %% full momentum step
            switch model
                case 'laplace'
                    dEda = dEdX_laplace(a(:), X, phi, sigma, lambda);
                case 'gauss'
                    dEda = dEdX_gauss(a(:), X, phi, sigma);
            end
            p = p - langevin_epsilon*dEda;

            %% final half step
            a = a + (1/2)*langevin_epsilon*p;
            p = -p;

            %% MH criteria
            switch model
                case 'laplace'
                    Hold = E_laplace(ainit, X, phi, sigma, lambda) + ...
                           1/2*sum(pinit.^2, 1);
                    Hnew = E_laplace(a, X, phi, sigma, lambda) + ...
                           1/2*sum(p.^2, 1);
                case 'gauss'
                    Hold = E_gauss(ainit, X, phi, sigma, lambda) + ...
                           1/2*sum(pinit.^2, 1);
                    Hnew = E_gauss(a, X, phi, sigma, lambda) + ...
                           1/2*sum(p.^2, 1);
            end
            pacc = exp( Hold - Hnew );
            rej = find(pacc < rand(size(pacc)));
            a(:,rej) = ainit(:,rej);
            p(:,rej) = pinit(:,rej);
            fprintf('\rt %d step %d reject frac %f', t, i, ...
                    length(rej) / size(a,2) );

            %% corrupt momentum
            R = randn(M, B);
            p = -sqrt(1-langevin_beta)*p + sqrt(langevin_beta)*R;
        end
    end

    fprintf('\n');


    time_inf = toc;

    p_log = [p_log, sqrt(sum(p(:).^2)/B/M)];
    switch model
        case 'laplace'
            energy = E_laplace(a(:), X, phi, sigma, lambda);
        case 'gauss'
            energy = E_gauss(a(:), X, phi, sigma);
    end
    energies_log = [energies_log, sum(energy)];

    if test_every == 1 || mod(update,test_every) == 0
        fig_ais

        loglike_ais_log = [loglike_ais_log loglike_ais];
    end


    %% m-step (parameter update)
    tic();

    E = X - phi*a;
    snr = 10 * log10(sum(X(:).^2) / sum(E(:).^2));

    if learn_sigma
        % update sigma
        sigma = E * E' / size(E,2);
    end
    
    %% use a newton step to optimize the quadratic energy in phi
    H = inv(a*a'/B);
    g = -E*a'/B;
    d = -g*H';
    phi = phi + d;

    len = sqrt(sum(sum(d.^2))/M);
    updatelength_log = [updatelength_log, len];

    time_lrn = toc();

    %% display

    if display_every == 1 || mod(update,display_every) == 0

        sfigure(6); clf;
        phi_mag = sqrt(sum(phi.^2, 1));
        [junk, ord] = sort(phi_mag);
        bar(phi_mag(ord));
        axis tight;
        title('basis lengths (sorted)');

        EI = V*phi*a(:,1);
        E = V*X(:,1) - EI;

        sc = max([abs(V * X(:,1)) ; abs(EI(:)) ; abs(E(:))]);

        sfigure(4);
        colormap(gray);
        subplot(1,3,1);
            imagesc(reshape(V*X(:,1), Jsz, Jsz), [-sc sc]);
            axis image; title('original');
        subplot(1,3,2);
            imagesc(reshape(EI, Jsz, Jsz), [-sc sc]);
            axis image; title('reconstructed');
        subplot(1,3,3);
            imagesc(reshape(E, Jsz, Jsz), [-sc sc]);
            axis image; title('difference');

        sfigure(5);
        bar(a(:,1)); axis tight;

        sfigure(1);
        colormap(gray);
        clim = max(max(abs(V*phi)));
        array_ord = render_network(V*phi(:, ord), Mrows, clim);
        imagesc(array_ord, [-1 1]);
        axis image off;

        if learn_sigma
            sfigure(70);
            imagesc(sigma); colorbar;
            axis image off;
            title('sigma (reconstruction covariance matrix)');
        end

        sfigure(2);
        subplot(2,1,1); hist(a(:), 100); title('a');
        subplot(2,1,2); hist(randn(prod(size(a)),1), 100); title('normal');



        %% learning param dynamics
        sfigure(60); clf;

        subplot(4,1,1);
            semilogy(updatelength_log, 'r-');
            axis tight;
            title('update step length');
            xlabel('learning step');
            grid on;

        subplot(4,1,2);
            plot(p_log, 'g-');
            axis tight;
            title('Langevin momentum length');
            grid on;

        subplot(4,1,3);
            plot(energies_log, 'b-');
            axis tight;
            title('energy history');
            grid on;

        if update >= test_every
            subplot(4,1,4);
            plot(test_every:test_every:update, loglike_ais_log, '.-', ...
                 test_every:test_every:update, mean(loglike_ais_log), '.--' );
            axis tight;
            title('loglike history');
            xlabel('learning step');
            grid on;
        end


        drawnow;
    end

    if save_every == 1 || mod(update,save_every) == 0
        [sucess,msg,msgid] = mkdir(sprintf('state/%s', paramstr));

        clim = max(max(abs(V*phi)));
        array = render_network(V*phi, Mrows, clim);
        array_frame = uint8(255*((array+1)/2)+1);

        imwrite(array_frame, ...
            sprintf('state/%s/bf_up=%06d.png',paramstr,update), 'png');
        eval(sprintf('save state/%s/phi.mat phi',paramstr));

        filename = sprintf('state/%s/bf_lengths.png',paramstr);
        saveas(6, filename);

    end

    fprintf('%s', paramstr);
    fprintf(' update %d snr %.4f', update, snr);
    fprintf(' len %.4f', len);
    fprintf(' a_var %.4f', var(a(:)) );
    fprintf(' a_abs %.4f', mean(abs(a(:))) );
    fprintf(' p_var %.4f', var(p(:)) );
    fprintf(' inf %.4fs lrn %.4fs', time_inf, time_lrn);
    fprintf('\n');

    update = update + 1;
end

eval(sprintf('save state/%s/matlab_up=%06d.mat', paramstr, update)); 


