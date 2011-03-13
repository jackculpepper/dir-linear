
switch datasource
    case 'vanHat'
        vHwidth = 1536;
        vHheight = 1024;
        
        [sucess,msg,msgid] = mkdir(sprintf('cache', paramstr));
        filename = sprintf('cache/Xtestr_vanHat_%dx%d.mat', J, Btest );

        if exist( filename, 'file')
            cache = load(filename);
            Xtestr = cache.Xtestr;
        else
            Xtestr = zeros(J,Btest);  
            file_list = importdata( '../data/vanHateren/image_list_test.txt');
            i_image = 0;
            for b = 1:Btest
                if mod(b-1, patches_per_image_test) == 0
                    % load a new image
                    imgfilename = file_list{floor(rand*length(file_list))+1};
                    imgfilename = ['../data/vanHateren/iml00001-04212/', imgfilename];
                    f1 = fopen(imgfilename,'rb','ieee-be');
                    img = fread(f1,[vHwidth,vHheight],'uint16');
                    fclose(f1);
                    img = img';
                    img = img / std(img(:));
                end

                r = buff + ceil((vHheight-Jsz-2*buff)*rand);
                c = buff + ceil((vHwidth-Jsz-2*buff)*rand);

                Xtestr(:,b) = reshape(img(r:r+Jsz-1,c:c+Jsz-1), J, 1);
                fprintf('\r%d / %d', b, Btest);
            end

            %% save to cache
            save(filename, 'Xtestr');
        end

        Xtestrl = Xtestr;
        if log_images
            Xtestrl = log( Xtestrl );
        end

        Xtestrl = bsxfun( @plus, Xtestrl, -mean( Xtestrl ) );

        % PCA project it onto L components
        Xtest = V'*Xtestrl;
        Xtest = diag(sqrt(1./diag(D))) * Xtest;
    case 'images'
        [sucess,msg,msgid] = mkdir(sprintf('cache', paramstr));
 
        filename = sprintf('cache/Xtest_images_%dx%dx%d.mat', J, L, Btest);
        if exist( filename, 'file')
            cache = load(filename);
            Xtest = cache.Xtest;
        else
            
            load ../data/IMAGES_RAW.mat
            [Nsz,Nsz,K] = size(IMAGESr);
                
            Xtestr = zeros(J,Btest);  

            % extract subimages at random
            fprintf('selecting image patches ..\n');
            for b = 1:Btest
                i = 10; %% select image 10 for testing
                r = buff + ceil((Nsz-Jsz-2*buff)*rand);
                c = buff + ceil((Nsz-Jsz-2*buff)*rand);

                Xtestr(:,b) = reshape(IMAGESr(r:r+Jsz-1,c:c+Jsz-1,i), J, 1);
        
                Xtestr(:,b) = Xtestr(:,b) - mean(Xtestr(:,b));

                fprintf('\r%d / %d', b, Btest);
            end
            fprintf('\n');

            Xtest = V'*Xtestr;
            Xtest = diag(sqrt(1./diag(D))) * Xtest;

            %% save to cache
            save(filename, 'Xtest');
        end

end

