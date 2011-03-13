switch datasource
    case 'vanHat'
        
        [sucess,msg,msgid] = mkdir(sprintf('cache', paramstr));
      
        filename = sprintf('cache/Xr_vanHat_%dx%d.mat', J, B);
        if exist( filename, 'file')
            cache = load(filename);
            Xr = cache.Xr;
        else
            vHwidth = 1536;
            vHheight = 1024;

            Xr = zeros(J,B);

            file_list = importdata( '../data/vanHateren/image_list_train.txt');
            fprintf('selecting image patches ..\n');

            [val,ord] = sort(randn(1,B));

            for b = 1:B
                if mod(b-1, patches_per_image) == 0
                    % load a new image
                    imgfilename = ['../data/vanHateren/iml00001-04212/', ...
                                   file_list{floor(rand*length(file_list))+1}];

                    f1 = fopen(imgfilename,'rb','ieee-be');
                    img = fread(f1,[vHwidth,vHheight],'uint16');
                    fclose(f1);

                    img = img';
                    img = img / std(img(:));

                    figure(1);
                    colormap(gray);
                    imagesc(img);
                    drawnow;
                end

                r = buff + ceil((vHheight-Jsz-2*buff)*rand);
                c = buff + ceil((vHwidth-Jsz-2*buff)*rand);

                Xr(:,ord(b)) = reshape(img(r:r+Jsz-1,c:c+Jsz-1), J, 1);
                fprintf('\r%d / %d', bb, B);
            end

            %% save to cache
            save(filename, 'Xr');
        end

        Xrl = Xr;
        if log_images ; Xrl = log( Xrl +1e-6 ); end

        %% subtract the mean
        Xrl = bsxfun(@plus, Xrl, -mean(Xrl));

        %% project onto L PCA components
        C = Xrl*Xrl' / B;
        [V,D] = eigs(C, L);
        X = V'*Xrl;

        %% whiten
        X = diag(sqrt(1./diag(D))) * X;

        clear Xr
        clear Xrl
    case 'images'
        [sucess,msg,msgid] = mkdir(sprintf('cache', paramstr));

        filename = sprintf('cache/X_images_%dx%dx%d.mat', J, L, B);
        if exist( filename, 'file')
            cache = load(filename);
            X = cache.X;
            V = cache.V;
            D = cache.D;
        else
            load ../data/IMAGES_RAW.mat
            [Nsz,Nsz,K] = size(IMAGESr);

            Xr = zeros(J,B);

            % extract subimages at random
            fprintf('selecting image patches ..\n');
            for b = 1:B
                i = ceil( (K-1)*rand );  %% leave image 10 for testing
                r = buff + ceil((Nsz-Jsz-2*buff)*rand);
                c = buff + ceil((Nsz-Jsz-2*buff)*rand);

                Xr(:,b) = reshape(IMAGESr(r:r+Jsz-1,c:c+Jsz-1,i), J, 1);

                Xr(:,b) = Xr(:,b) - mean(Xr(:,b));

                fprintf('\r%d / %d', b, B);
            end
            fprintf('\n');

            C = Xr*Xr' / B;
            [V,D] = eigs(C, L);
            X = V'*Xr;
            X = diag(sqrt(1./diag(D))) * X;

            %% save to cache
            save(filename, 'X', 'V', 'D');

        end


end


