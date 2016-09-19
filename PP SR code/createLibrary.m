function libraryPatch = createLibrary(img, patchSize, K)

[h, w] = size(img);
N = h*w;
side = (patchSize - 1)/2;
img = padarray(img, [side, side], 'replicate');

if K > N
    disp('Please make sure the number of library patches requested is lesser than number pixels in your image.');
else
    libraryPatch = zeros(patchSize, patchSize, K);
    randIndices = randperm(N, K);
    
    for k = 1:K
        i = floor((randIndices(k)-1)/w) + 1;
        j = randIndices(k) - (i-1)*w;
        libraryPatch(:, :, k) = img(i:i+(2*side), j:j+(2*side));
    end
end

