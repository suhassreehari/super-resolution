function high_res_image = cubicInterpolate(low_res_image, factor)

[h, w] = size(low_res_image);
high_res_image = zeros(h, w);
temp = griddedInterpolant(double(low_res_image), 'cubic');

for i = 1:h
    for j = 1:w
        for k = 1:factor
            for l = 1:factor
                high_res_image((i-1)*factor + k, (j-1)*factor + l) = temp(i + (k-1)/factor, j + (l-1)/factor);
            end
        end
    end
end
