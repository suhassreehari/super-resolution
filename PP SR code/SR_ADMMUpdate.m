function [map_image] = SR_ADMMUpdate(map_image, Mask, params, mu, R)

%Function that performs random order ICD with a quadratic penalty
%on the input (for Augmented lagrangian)
%Inputs map_image : initial value of image to be reconstructed
%           Mask  : Binary Mask of coordinate where measurements are done
%           params: hold the parameter for Augmented Lagrangian variables
% Outputs  map_image : Final reconstruction


TempK = params.v - params.u;

[m, n] = size(map_image);

for row = 1:m
    for col = 1:n
        CErow = (ceil(row/R)-1)*R + 1;
        CEcol = (ceil(col/R)-1)*R + 1;
        if(Mask(row, col) == 1)
            update = map_image(row, col);
        else
            update = TempK(row, col) + map_image(CErow, CEcol) - mu(CErow, CEcol);
        end
        
        if(update > 0)
            map_image(row, col) = update;
        else
            map_image(row, col) = 0;
        end
    end
end