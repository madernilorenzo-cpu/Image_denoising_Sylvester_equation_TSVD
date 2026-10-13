function img = load_example_image(n, m)
%LOAD_EXAMPLE_IMAGE Loads a test image of size N x M in double [0,1].
%
%   IMG = LOAD_EXAMPLE_IMAGE(N, M) looks for a standard MATLAB image
%   ('peppers.png'); if unavailable it generates a synthetic image
%   (blurred Gaussian blobs). No Image Processing Toolbox required:
%   grayscale conversion is done manually and no imresize is used, only
%   cropping.

try
    raw = imread('peppers.png');
    img = im2double01(raw);
    img = img(1:min(n, size(img,1)), 1:min(m, size(img,2)));
    img = padTo(img, n, m);
catch
    img = synthetic_image(n, m);
end
end

function x = im2double01(x)
%IM2DOUBLE01 Conversion to double [0,1] (uint8/uint16/double).
if isa(x, 'uint8')
    x = double(x) / 255;
elseif isa(x, 'uint16')
    x = double(x) / 65535;
end
if size(x, 3) == 3
    % ITU-R BT.601 conversion weights (equivalent to rgb2gray)
    x = 0.2989*x(:,:,1) + 0.5870*x(:,:,2) + 0.1140*x(:,:,3);
end
end

function x = padTo(x, n, m)
%PADTO Brings the image to size n x m by replicating the borders.
[nn, mm] = size(x);
if nn < n, x = padarrayReplicate(x, n - nn, 1); end
if mm < m, x = padarrayReplicate(x, m - mm, 2); end
end

function x = padarrayReplicate(x, k, dim)
%PADARRAYREPLICATE Replicates the border k times along dimension dim.
if dim == 1
    x = [x; repmat(x(end, :), k, 1)];
else
    x = [x, repmat(x(:, end), 1, k)];
end
end

function img = synthetic_image(n, m)
%SYNTHETIC_IMAGE Smooth synthetic image, values in [0,1].
[X, Y] = meshgrid(linspace(-1, 1, m), linspace(-1, 1, n));
img = 0.4 + 0.25*sin(2*X) .* cos(1.5*Y);
centers = [-0.5 -0.4; 0.4 0.3; 0.1 -0.6; -0.3 0.5];
amps = [0.5 0.4 0.3 0.35];
for i = 1:size(centers, 1)
    img = img + amps(i) * exp(-((X-centers(i,2)).^2 + (Y-centers(i,1)).^2) / 0.08);
end
img = max(0, min(1, img));
end
