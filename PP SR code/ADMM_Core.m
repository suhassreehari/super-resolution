function [map_image, params, eps_primal, eps_dual] = ADMM_Core(map_image, sino, geom, Amatrix, params, kparams, mode, forwardModel, libraryPatch, patchSize, R, HR_img)

%Performs ADMM based minimization of min_{x,v} (1/2)||y-Ax||^2 + s(v)
% Inputs
% map_image is an inital estimate for x
% sino is a structure which contains all the sinogram information
% geom contains the geometry information like number of pixels etc.
% Amatrix = A stored as a sparse matrix
% params contains parameters associated with the optimization
% kparams contains prior model parameters
% mode represents what prior is used in s(v)
%    0 - kSVD
%    1 - BM3D
%    2 - TV
%    3 - PLOW
%    4 - qGGMRF
%    5 - Otsu's segmentation
%    6 - SMAP segmentation %DOES NOT WORK 5/4/2013
%    7 - MAP segmentation

iter = 1;

stop_crit = ones(1, params.max_iter);
eps_primal = zeros(1, params.max_iter);
eps_dual = zeros(1, params.max_iter);

[m, n] = size(map_image);

%% DEBUG STATEMENTS
if (params.verbose)
    RMSE(iter) = sqrt(sum((map_image(:) - params.original(:)).^2)/numel(params.original));
end


%% END OF DEBUG STATEMENTS

while (iter < params.max_iter && stop_crit(iter) > params.threshold)
    
    
    if(params.verbose)
        display(iter);
    end
    
    if(iter > 1) %After the first outer iteration just do 1 ICD iteration
        params.num_iter = 1;
        
        if(isfield(kparams,'initdict') && ischar(kparams.initdict))
            kparams.initdict=dict; %for kSVD use the previous dictionary to start
        end
    end
    
    
    switch(forwardModel)
        case 0
            map_image = DataTermOpt(map_image, sino,geom, params, Amatrix);
        case 1
            map_image = DataTermOptDeblur(map_image, sino, params, Amatrix);
        case 2
            map_image = DataTermOptInPaint(map_image, sino, params, Amatrix);
        case 3
            map_image = DataTermOpt3D(map_image, sino, geom, params, Amatrix);
        case 4
            [m, n] = size(map_image);
            mu = zeros(m, n);
            
            for row = 1:R:m
                for col = 1:R:n
                    mu(row, col) = mean(mean(map_image(row:row+R-1, col:col+R-1)));
                end
            end
            
            map_image = DataTermOptSR(map_image, sino, params, Amatrix, mu, R); % Super resolution
    end
    
    kparams.x = map_image + params.u;
    
    switch(mode)
        
        %%ksvd
        case 0
            kparams.x = (kparams.x./kparams.maxval).*255;%Normalize to 255
            [imout,dict] = ksvddenoise(kparams, 0);
            imout = (imout./255).*kparams.maxval;%Normalize this back to the regular scale
        case 1
            %%BM3D ; maxval is the max of noisy images (approximately)
            [~, imout] = BM3D(1, kparams.x/kparams.maxval, kparams.sigma);
            imout = imout.*kparams.maxval;%FOR BM3D this is needed
        case 2
            %TV denoising
            [imout, ~, ~, ~] = perform_tv_denoising(kparams.x, kparams);
            %imout=imout.*255;
        case 3
            kparams.x = (kparams.x./kparams.maxval).*255;%Normalize to 255
            Output = plowMatlab(zeros(m,n), kparams.x, kparams.sigma);
            imout = Output(:, :, end);
            imout = (imout./255).*kparams.maxval;%Normalize this back to the regular scale
        case 4
            imout = qGGMRFdenoise(kparams.x, kparams);
        case 5
            H = hamming(3)*hamming(3)';
            H = H./sum(H(:));
            kparams.x = imfilter(kparams.x,H);
            imout = OtsuSegmentationWrapper(kparams.x, kparams.num_class);
        case 6
            H = hamming(3)*hamming(3)';
            H = H./sum(H(:));
            kparams.x = imfilter(kparams.x,H);
            Temp = OtsuSegmentationWrapper(kparams.x, kparams.num_class);
            imout = SMAPSegmentationWrapper(Temp, kparams.num_class);
        case 7
            H = hamming(3)*hamming(3)';
            H = H./sum(H(:));
            Initial = imfilter(kparams.x,H);
            [Temp,~] = otsu(Initial,kparams.num_class);%Get a initial label set from Otsu
            [imout,~,~] = MAP_segmentation(kparams.x,kparams,Temp);
        case 8
            [~, imout] = bm4d(1,kparams.x/kparams.maxval,kparams.sigma);
            imout = imout.*kparams.maxval;%FOR BM4D this is needed
        case 9 % 2D NLM
            % kparams.x = (kparams.x./kparams.maxval).*255; % Normalize to 255
            if iter < 15
                I = double(kparams.x);
                %I = I/kparams.maxval;
                %imout = FastNonLocalMeans3D(I, kparams.sigma/255, params.beta);
                %imout = imout*kparams.maxval;
                
                %imout = NLM_ORACLE(I, HR_img, 0, 15+iter, 3, iStart, iEnd, jStart, jEnd);
                
                imout = NLM(I, 0, 12+iter, 5, 40);
            else
                imout = kparams.x;
            end
        case 10 % 2D NLM_ORACLE
            if iter < 10
                I = double(kparams.x);
                %imout = NLM_img2(I, 0, 14+iter, patchSize, 4, HR_img);
                [imout, weight] = NLM_ORACLE(I, libraryPatch, 0, 20+iter, patchSize);
            else
                I = double(kparams.x);
                %imout = I;
                [imout, weight] = NLM_ORACLE(I, libraryPatch, 0, 20+iter, patchSize, weight);
            end
        case 11 % p_LB-NLM (considers only p best patches)
            I = double(kparams.x);
            [imout, weight] = LBNLM_p(I, libraryPatch, 0, 12, patchSize);
    end
    
    prev_v = params.v;
    params.v = imout;
    
    params.u = params.u + (map_image - params.v); % update of u
    
    %eps_primal(iter) = sum(abs(params.v(:)-map_image(:)))./length(map_image(:)); %sum(abs(map_image(:)));
    
    %eps_dual(iter) = sum(abs(params.v(:)-prev_v(:)))./length(map_image(:)); %./sum((abs(prev_v(:))));
    
    eps_primal(iter) = sqrt(sum(params.v(:)-map_image(:)).^2)./length(map_image(:)); %sum(abs(map_image(:)));
    eps_dual(iter) = sqrt(sum(params.v(:)-prev_v(:)).^2)./sqrt(sum(params.u(:).^2)); %./sum((abs(prev_v(:))));
    
    stop_crit(iter) = (eps_primal(iter) + eps_dual(iter))/2;
    
    %% Code to vary penalty parmeter
    if(isfield(params,'tau') && isfield(params,'mu'))
        
        if (eps_primal(iter) > params.mu*eps_dual(iter))
            params.lambda= params.lambda*params.tau;
            params.u=params.u/params.tau;
            display('Adjusting scaling');
        else
            if (eps_dual(iter) > params.mu*eps_primal(iter))
                params.lambda= params.lambda/params.tau;
                params.u=params.u*params.tau;
                display('Adjusting scaling');
            end
        end
    end
    %% End of varying penalty parameter segment
    
    iter = iter + 1;
    
    if (params.verbose == 1)
        display('Average stopping criteria');
        display(stop_crit(iter-1));
        display('RMSE');
        display(RMSE(iter-1));
        
        if (mode == 4 && kparams.verbose == 1)
            %Debug code - cost computation for qGGMRF
            cost(iter) = 0;
            [sinonew] = forward_project_v2(map_image, sino, Amatrix);
            cost(iter) = (1/2)*sum(D(:).*(sino.counts(:) - sinonew(:)).^2);
            cost(iter) = cost(iter)+(params.beta*kparams.p*kparams.sigmax^kparams.p)*qGGMRFprior2DCost(map_image,kparams);
            
            if (cost(iter) - cost(iter-1)>0)
                display('Cost just increased!');
            end
        end
    end
    
    RMSE(iter) = sqrt(sum((map_image(:)-params.original(:)).^2)/numel(params.original));
end

params.RMSE = RMSE;
params.eps_primal = eps_primal;
params.eps_dual = eps_dual;

params.verbose = 0;

if (params.verbose == 1)
    if (mode == 0) %kSVD display final dictionary
        dictimg = showdict(dict,[1 1]*kparams.blocksize,round(sqrt(kparams.dictsize)),round(sqrt(kparams.dictsize)),'lines','highcontrast');
        figure; imshow(imresize(dictimg,2,'nearest'));
        title('Trained dictionary');
        display(iter);
        params.dictimg = dictimg;
    end
    if (mode == 4 && kparams.verbose == 1) %qGGMRF debug
        figure;
        plot(cost);
        cost_admm = cost;
        RMSE_ADMM = RMSE;
        save('ADMMCostv2.mat', 'cost_admm', 'RMSE');
    end
    
    %Plot the primal and dual convergence checks
    
    figure; plot(eps_primal(1:iter-1));% format_plot_for_publication(gcf);
    
    %title('||x-v||/||x||');
    %figure;
    %hold on;
    %plot(eps_dual(1:iter-1),'r'); % format_plot_for_publication(gcf);
    %h=legend('$\frac{1}{N}||x^{k}-v^{k}||$','$\frac{1}{N}||v^{k+1}-v^{k}||$');
    %set(h,'Interpreter','latex');
    xlabel('Iteration number (k)');
    %ylabel('RMS difference between $\hat{x}$ and $\hat{v}$');
    ylabel('Scaled primal residue');
    %title('||v^{k+1}-v^{k}||/||v^{k}||');
    
    figure; plot(eps_dual(1:iter-1));% format_plot_for_publication(gcf);
    xlabel('Iteration number (k)'); ylabel('Scaled dual residue');
end

