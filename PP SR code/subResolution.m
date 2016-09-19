function subResImg = subResolution(img, factor)
% This function artifically reduces resolution of an image by the
% user-defined factor

if nargin == 1
    factor = 4;
end

[h, w] = size(img);

subResImg = zeros(h/factor, w/factor);

for i = 1:factor:h
    for j = 1:factor:w
        subResImg((i-1)/factor + 1, (j-1)/factor + 1) = mean(mean(img(i:i+factor-1, j:j+factor-1)));
    end
end



end

