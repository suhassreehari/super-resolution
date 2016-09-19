function [Amatrix] = GenSamplingMask(image, indicator)
%Function that extracts a sampling mask from a image data set
%Inputs     image : An image from which certain pixels are "missing"
%           indicator  : The value used in the unmeasured locations 
% Outputs  A_matrix : A array containing locations which are measured

[I, J] = find(image ~= indicator);
Amatrix = [I, J];