function [Object, sino, geom, Amatrix] = data_creation(image_name, mask, sampled, InitialEst)


% InPainting data creation (Code by Suhas Sreehari, May 2015)
% -------------------------------------------------------------------------

% Input and parameters:

% Object = double(imread(image_name)); % [0, 255] scale

Object = InitialEst; % (since there is no ground truth available)

% percentage_fill = sampling; % fill percentage; default: 40
%sino.var_w = 100; % noise variance -- on a [0, 255] scale; default: 100
sino.var_w = 0;

% -------------------------------------------------------------------------

%sino.counts = imnoise(Object/255, 'gaussian', 0, sino.var_w/(255*255));
%sino.counts = sino.counts*255;
sino.counts = sampled;
[h, w] = size(Object);
geom.n_x = h;
geom.n_y = w;

% p_f = percentage_fill/100;

% mask = ones(h, w); % initilaizing to 1 (all pixels are by default measured)
% count = 0;

% for i = 1:h
%     for j = 1:w
%         if (rand > p_f) % pixel must not be measured
%             count = count + 1;
%             mask(i, j) = 0;
%             %Amatrix(count, 1) = i;
%             %Amatrix(count, 2) = j;
%             sino.counts(i, j) = 0;
%         end
%     end
% end


Amatrix = GenSamplingMask(mask, 0);

clearvars -except Object geom sino Amatrix mask

sampled_image = sino.counts;
