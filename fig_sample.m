

%% draw samples from the prior
Bgen = 100;

switch model
    case 'gauss'
        a_gen = randn(M, B);
    case 'laplace'
        a_gen = rand(M, B);
        a_gen = expinv(a_gen, 1/lambda);
end

%% permutation to select random data items
perm = randperm(B);

%% pass through generative function
a_gen_mod = a_gen(:,1:Bgen);
a_gen_pst = a_gen(:,perm(1:Bgen));

%% pass through PCA basis
Xhat_mod = V * phi * a_gen_mod;
Xhat_act = V * X(:,perm(1:Bgen));
Xhat_pst = V * phi * a_gen_pst;

clim_mod = max(abs(Xhat_mod(:)));
clim_act = max(abs(Xhat_act(:)));
clim_pst = max(abs(Xhat_pst(:)));

array_Xhat_mod = render_network(Xhat_mod, sqrt(Bgen), clim_mod);
array_Xhat_act = render_network(Xhat_act, sqrt(Bgen), clim_act);
array_Xhat_pst = render_network(Xhat_pst, sqrt(Bgen), clim_pst);


sfigure(124); colormap(gray);
imagesc(array_Xhat_mod, [-1 1]); axis image off;
title(sprintf('Draws from full model, scale = %.4f', clim_mod));

sfigure(126); colormap(gray);
imagesc(array_Xhat_act, [-1 1]); axis image off;
title(sprintf('Actual image patches, scale = %.4f', clim_act));

sfigure(127); colormap(gray);
imagesc(array_Xhat_pst, [-1 1]); axis image off;
title(sprintf('Draws from the posterior, scale = %.4f', clim_pst));

array_Xhat_mod_frame = uint8(255*((array_Xhat_mod+1)/2)+1);
array_Xhat_act_frame = uint8(255*((array_Xhat_act+1)/2)+1);
array_Xhat_pst_frame = uint8(255*((array_Xhat_pst+1)/2)+1);

array_Xhat_mod_frame = imresize(array_Xhat_mod_frame, 8, 'nearest');
array_Xhat_act_frame = imresize(array_Xhat_act_frame, 8, 'nearest');
array_Xhat_pst_frame = imresize(array_Xhat_pst_frame, 8, 'nearest');

if exist(sprintf('state/%s', paramstr), 'dir')
    imwrite(array_Xhat_mod_frame, sprintf('state/%s/fig_Xhat_mod.png', ...
        paramstr), 'png');
    imwrite(array_Xhat_act_frame, sprintf('state/%s/fig_Xhat_act.png', ...
        paramstr), 'png');
    imwrite(array_Xhat_pst_frame, sprintf('state/%s/fig_Xhat_pst.png', ...
        paramstr), 'png');
end


