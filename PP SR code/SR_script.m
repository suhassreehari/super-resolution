
N = 1; % number of images
patchSize = 15; % patch size for NLM-based methods; default: 5

image_name = cell(N, 1);
%image_name{1} = 'Images/Gold_atoms_4xsubRes.tif';
image_name{1} = 'Images/E_arial/E_LR6x.png';
library_image_name{1} = 'Images/E_arial/lib_HR.png';
% ground_truth{1} = 'Images/Crystal_L1.tif';

output_folder_name = cell(N, 1);
output_folder_name{1} = 'Results/E_arial';

method = [1, 10, 0];
method_iter = [50, 2, 20];
method_B = [1, 1, 7];

R = 6;

initializationMethod = 1;  % 0: Shepard, 1: Cubic, 2: use your own initialization


for i = 1:N
    img = double(imread(image_name{i}));
    [h, w] = size(img);
    
    NL1 = 0;
    NL2 = 3000; % Default: 5000
    HR_img = double(imread(library_image_name{i})); % High-Resolution image
    
    library1 = createLibrary(img, patchSize, NL1);
    library2 = createLibrary(HR_img, patchSize, NL2);
    
    library = zeros(patchSize, patchSize, NL1+NL2);
    library(:, :, 1:NL1) = library1(:, :, :);
    library(:, :, NL1+1:NL1+NL2) = library2(:, :, :);
    
    
    % 4x super-resolution case
    if R == 4
        
        mask = zeros(R*h, R*w);
        mask(1:R:R*h, 1:R:R*w) = 1;
        sampled = zeros(R*h, R*w);
        sampled(1:R:R*h, 1:R:R*w) = img(1:h, 1:w);
        
        if (initializationMethod == 0) % Shepard's interpolation
            InitialEst = shepard_initialize(sampled, mask, 21, 2);
            InitialEst = InitialEst(1:R*h, 1:R*w);
            InitialEst(1:R:R*h, 1:R:R*w) = img(1:h, 1:w);
            
        elseif (initializationMethod == 1) % Cubic interpolation
            InitialEst = cubicInterpolate(img, R);
        elseif (initializationMethod == 2) % Load your own initialization
            InitialEst = double(imread('Results/a3/4_NLM_ORACLE_size11.png'));
        end
        
        % 6x super-resolution case
    elseif R == 6
        
        mask = zeros(R*h, R*w);
        mask(1:R:R*h, 1:R:R*w) = 1;
        sampled = zeros(R*h, R*w);
        sampled(1:R:R*h, 1:R:R*w) = img(1:h, 1:w);
        
        if (initializationMethod == 0) % Shepard's interpolation
            InitialEst = shepard_initialize(sampled, mask, 21, 2);
            InitialEst = InitialEst(1:R*h, 1:R*w);
            InitialEst(1:R:R*h, 1:R:R*w) = img(1:h, 1:w);
            
        elseif (initializationMethod == 1) % Cubic interpolation
            InitialEst = cubicInterpolate(img, R);
        elseif (initializationMethod == 2) % Load your own initialization
            InitialEst = double(imread('Results/a3/4_NLM_ORACLE_size11.png'));
        end
        
    end
    
    say = sprintf('Starting image # %d/%d...', i, N);
    disp(say);
    
    for k = 2:2
        close all
        clc
        
        [map_image, params, primal_residue, dual_residue] = PlugAndPlayPriorFunction(image_name{i}, output_folder_name{i}, R, method(k), method_iter(k), method_B(k), InitialEst, mask, sampled, library, patchSize, R, HR_img);
        
        say = sprintf('Done with combination # %d/%d', k, K);
        disp(say);
    end
    
    say = sprintf('Done with image # %d/%d', i, N);
    disp(say);
end