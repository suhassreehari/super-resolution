function [map_image, params, primal_residue, dual_residue] = PlugAndPlayPriorFunction(image_name, output_folder_name, super_res, method, method_iter, method_B, InitialEst, mask, sampled, library, patchSize, R, HR_img)

% Plug and Play priors script file
% To solve (1/2)||y-Ax||_{D}^2+ beta*prior(x)
% sino - (Tomographic projection) data ntheta X nr
% geom - object dimensions
% lambda - lagrange multiplier
% beta - regularization 
% u - augmented lagrangian vector nx X ny (Object size)
% v - temporary image used in AL method nx X ny

if nargin == 0
    disp('Error: need at least the name of the image.');
elseif nargin == 1
    super_res = 4;
    method = 1; % BM3D by default
elseif nargin == 2
    method = 1;
end

%load('CameramanData40%StdDev10.mat');

[Object, sino, geom, Amatrix] = data_creation(image_name, mask, sampled, InitialEst); % line added by Suhas Sreehari.

out_file_name = sprintf('/%d_sampling.png', super_res);
out_file_name = strcat(output_folder_name, out_file_name);
imwrite(uint8(sino.counts), out_file_name);

out_file_name = sprintf('/%d_shepard.png', super_res);
out_file_name = strcat(output_folder_name, out_file_name);
imwrite(uint8(InitialEst), out_file_name);

out_file_name = sprintf('/%d_mask.png', super_res);
out_file_name = strcat(output_folder_name, out_file_name);
imwrite(uint8(mask*255), out_file_name);

WriteFile = 0;
Display = 0;
ScaleMin = 0;
ScaleMax = 255; % Dynamic range of input data
NumRes = 1; % This is for multi-resolution - ONLY CHANGE FOR TOMOGRPHY

if(WriteFile == 1)
    Path = '/Users/suhas/Documents/';
    Filename = 'DR.eps';
    FullPath = strcat(Path, Filename);
end


for res = 1:NumRes
    
    new_geom = geom;
    new_geom.n_x = new_geom.n_x/(2^(NumRes-res));
    new_geom.n_y = new_geom.n_y/(2^(NumRes-res));
    
    m = new_geom.n_y;
    n = new_geom.n_x;
    
    ObjectSize = max(m, n);
    ObjectRes = 1; % Size of object voxels relative to detector pixel
    DetDim = 1; % 1.5*1000;%um - micron
    
    %% Prior Model Selection
    
    %Choose the denoising algorithm 0 - kSVD, 1-BM3D, 2-TV, 3 - PLOW, 4 -
    %qGGMRF, 5 - Otsu's segmentation, 6 - SMAP segmentation - BUGGY!! TODO: Repair
    %7 - MAP segmentation, 8 - BM4D, 9 - 2D NLM, 10 - 2D library NLM.
    mode = method;
    
    %Chooise forward model 0 - 2-D tomography 1 - deblurring 2 -
    %inpainting/interpolation
    %3 - 3D tomography, 4 - 2D super resolution.
    forwardModel = 4;
    
    params.max_iter = method_iter;
    B = method_B;
    
    params.threshold = 1e-3;% Stopping criteria for algorithm
    
    %% Other paramters - ICD
    
    %%%%%%%%%%% - USER DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
    params.lambda = 1/20; % Default 1/20
    params.beta = B*params.lambda; % default: 30*params.lambda
    
    %%%%%%%%%%% END OF USER DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
    
    params.num_iter = 1; % Number of ICD inner loops per ADMM iter
    params.filter = [1/12, 1/6, 1/12; 1/6, 0, 1/6; 1/12, 1/6, 1/12];
    params.u = zeros(m, n); % Augmented Lagrange vector
    params.v = zeros(m, n); % Auxilary variable
    params.verbose = 1;
    params.original = imresize(Object, [m, n]); % Phantom used
    
    params.xray = 0;
    
    %% Denoising parameter (std. dev) for ALL algorithms
    kparams.sigma = sqrt(params.beta/params.lambda);
    
    %% Denoising parameters - KSVD and BM3D
    if (mode == 0 || mode == 1 || mode == 9)
        %kparams.x = (map_image./max(map_image(:))).*255;
        kparams.blocksize = 2; %4;
        kparams.dictsize = 512; %256
        kparams.maxval= ScaleMax;
        kparams.trainnum = 3600; %40000;%40000;%%40000;
        kparams.iternum = 10;
        kparams.memusage = 'high';
    end
    %% These are only for TV denoising
    if (mode == 2)
        kparams.verb = 0;
        kparams.display = 0;
        kparams.niter = 100; % number of iterations
        kparams.c_TV = 0.0998; %0.0998;%10;%0.0998;%0.0998;%1;%.0998;
        %kparams.etgt = (kparams.sigma/255)*(max(m,n));
        kparams.lambda = kparams.c_TV*kparams.sigma^2; %1e-4;%1e-4;%kparams.c_TV*kparams.sigma^2;%2./kparams.sigma^2; % initial regularization
    end
    
    %% PLOW
    if (mode == 3)
        kparams.maxval = ScaleMax;
    end
    
    %% qGGMRF
    if (mode == 4)
        kparams.p = 2;
        kparams.q = 1.2;
        kparams.c = 1e-2; % 1e-3;
        kparams.sigmax = 0.4; %0.29; % 0.00035*2;%0.29*256;%0.29;%006/2;%.00035;%%%%.1*.00035;%.%.0012/2;%.0002;;%0.5940*(params.lambda^(1/kparams.p));
        kparams.niter = 10; % CHANGE THIS TO AVOID CONFLICT
        kparams.verbose = 0;
        display('Bad variable names for num iter - denoising!');
    end
    %% General Segmentation parameter
    if (mode == 7)
        kparams.num_class = 6; %3
        
        %% MAP Segmentation
        kparams.max_iter = 10;
        kparams.filter = params.filter;
        kparams.beta = 4; % 3000;
        kparams.debug = 0;
        kparams.sigma_sq = (1/(params.lambda));
        kparams.rand_ord = 0; % Regular ICM or random order ICM
    end
    
    if (mode == 8)
        kparams.maxval = ScaleMax;
    end
    
    if (mode == 9) || (mode == 10)
        kparams.maxval = ScaleMax;
    end
        
    %% Augmented Lagrangian Iterations
    
    
    start = tic;
    
    %initialization for inpainting
    mask = zeros(size(sino.counts));
    for i = 1:length(Amatrix)
        row = Amatrix(i,1);
        col=Amatrix(i,2); mask(row,col) = 1;
    end
    
    params.v = InitialEst;
    [map_image, params, primal_residue, dual_residue] = ADMM_Core(InitialEst, sino, new_geom, Amatrix, params, kparams, mode, forwardModel, library, patchSize, R, HR_img);
    
    elapsed = toc(start);
    
    
end

display(elapsed);

RMSE = sqrt(sum(sum((map_image - Object).^2))/numel(Object));

figure; imagesc(map_image, [ScaleMin, ScaleMax]); axis image; colormap(gray);
% format_image_for_publication(gcf);

if (Display == 1)
    
    figure;
    imagesc(InitialEst, [ScaleMin, ScaleMax]); axis image; colormap(gray)
    colorbar('Eastoutside');
    format_image_for_publication(gcf);
    title('Initial reconstruction');
    
    figure;
    imagesc(map_image, [ScaleMin, ScaleMax]); axis image; colormap(gray)
    colorbar('Eastoutside');
    format_image_for_publication(gcf);
    title(strcat('RMSE = ', num2str(RMSE), ' \beta = ', num2str(params.beta)));
    
    figure;
    imagesc(params.v, [ScaleMin, ScaleMax]); axis image; colormap(gray)
    colorbar('Eastoutside');
    format_image_for_publication(gcf);
    title('v');
    
    if(params.verbose)
        figure;
        plot(params.RMSE); format_plot_for_publication(gcf);
        xlabel('Iteration number'); ylabel('RMSE');
    end
    
end

%% Write file to TechReport path

if(WriteFile == 1)
    figure;
    imshow(uint8(map_image));
    format_image_for_publication(gcf);
    export_fig(FullPath, '-eps');
end

% m = medfilt2(map_image);

clear out_file_name;

% if method == 0 && mod(sampling, 5) == 0
%     out_file_name = sprintf('/%d_NLM_rmse_%d.png', sampling, params.RMSE(1));
% elseif method == 0 && mod(sampling, 5) ~= 0
%     out_file_name = sprintf('/%d_DSG-NLM_rmse_%d.png', sampling - mod(sampling, 5), params.RMSE(1));
% elseif method == 1
%     out_file_name = sprintf('/%d_BM3D_rmse_%d.png', sampling, params.RMSE(1));
% end

% if mod(method_k - 1, 3) == 0
%     out_file_name = sprintf('/%d_BM3D_rmse_%d.png', sampling, params.RMSE(1));
% elseif mod(method_k - 2, 3) == 0
%     out_file_name = sprintf('/%d_DSG-NLM_rmse_%d.png', sampling - mod(sampling, 5), params.RMSE(1));
% elseif mod(method_k, 3) == 0
%     out_file_name = sprintf('/%d_NLM_rmse_%d.png', sampling, params.RMSE(1));
% end

if method == 0
    out_file_name = sprintf('/%d_KSVD_rmse_%d.png', super_res, params.RMSE(1));
elseif method == 1
    out_file_name = sprintf('/%d_BM3D_rmse_%d.png', super_res, params.RMSE(1));
elseif method == 4
    out_file_name = sprintf('/%d_qGGMRF_rmse_%d.png', super_res, params.RMSE(1));
elseif method == 8
    out_file_name = sprintf('/%d_NLM_rmse_%d.png', super_res, params.RMSE(1));
elseif method == 9
    out_file_name = sprintf('/%d_DSG-NLM_rmse_%d.png', super_res, params.RMSE(1));
elseif method == 10
    out_file_name = sprintf('/%d_NLM_ORACLE_rmse_%d.png', super_res, params.RMSE(1));
%     out_file_name2 = sprintf('/%d.png', 1);
elseif method == 11
    out_file_name = sprintf('/%d_LBNLM_p_rmse_%d.png', super_res, params.RMSE(1));
end

out_file_name = strcat(output_folder_name, out_file_name);
% imwrite(uint8(map_image), out_file_name);
imwrite(uint8(params.v), out_file_name);

% out_file_name2 = strcat(output_folder_name, out_file_name2);
% imwrite(uint8(params.v), out_file_name2);

