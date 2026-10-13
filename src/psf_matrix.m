function P = psf_matrix(n, type)
%PSF_MATRIX Point spread function matrix (Section 3.3.1).
%
%   P = PSF_MATRIX(N, TYPE) builds the N x N matrix associated with the
%   separable, space-invariant PSF P = P_y (kron) P_x.
%   The available choices are those listed in Section 3.3.1:
%
%     'identity' : P = I                              (choice 1)
%     'average'  : band blur, r = 4                   (choice 2)
%     'gaussian' : Gaussian band, r = 35, sigma = 5   (choice 3)
%
%   Columns are normalized to sum to one.

switch lower(type)
    case 'identity'
        P = eye(n);
        return
    case 'average'
        r = 4;
        e = (-(r-1):(r-1)).';
        w = ones(size(e)) / (2*r-1);
    case 'gaussian'
        r = 35;
        s = 5;
        e = (-r:r).';
        w = 0.5 * exp(-(e.^2) / (2*s^2));
    otherwise
        error('psf_matrix:unknownType', 'Unknown PSF type: %s', type);
end

P = zeros(n);
idx = (1:n).';
for k = 1:numel(e)
    d = e(k);
    valid = (idx + d >= 1) & (idx + d <= n);
    P(idx(valid) + d, idx(valid)) = w(k);
end
P = P ./ max(sum(P, 1), eps);
end
