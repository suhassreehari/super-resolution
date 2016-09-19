function [new_image] = DataTermOptSR(map_image, sino, params, List, mu, R)

%Optimizes the function : data mismatch term + augmented vector
% min_x ((1/2)||y-Ax||_D^{2} + lambda ||x-(v-rho)||^2 )
%Inputs: map_image : Initial value of x 
%        sino : sinogram structure containing data with missing pixels
%        geom : geometry structure contianing geom info of object 
%        params : prior model parameters + algorithm parameters structure
%        List : List of pixels where measurement is available 


%Compute initial error image

%F = fspecial('gaussian',5);

Mask = zeros(size(sino.counts));
[ListLen, ~] = size(List);
for i = 1:ListLen
    row = List(i, 1);
    col = List(i, 2);
    Mask(row,col) = 1;
end

for k = 1:params.num_iter          
    % Update all the pixels once
    [map_image] = SR_ADMMUpdate(map_image, Mask, params, mu, R);   
end

new_image = map_image;
