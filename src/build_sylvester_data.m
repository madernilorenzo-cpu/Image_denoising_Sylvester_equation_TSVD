function [A, B, C, mu, sigma_f2] = build_sylvester_data(G, opts)
%BUILD_SYLVESTER_DATA Builds A, B, C for the denoising problem.
%
%   [A, B, C, MU, SIGMA_F2] = BUILD_SYLVESTER_DATA(G, OPTS) builds the
%   terms of the Sylvester equation
%
%       A * Fhat + Fhat * B = C
%
%   associated with the denoising problem described in Sections 2.2 and
%   3.3.1 of the thesis, starting from the degraded image G (N x M).
%
%   Model (eq. (2.15)-(2.16), with Sigma_z = SIGMA_F2 * R_z):
%
%       A = SIGMA2 * (P_y * Sigma_y)^{-1} = (SIGMA2/SIGMA_F2) * (P_y*R_y)\I
%       B = P_x * Sigma_x               = SIGMA_F2 * (P_x * R_x)
%       C = P_y^{-1} * G * (P_x Sigma_x)
%
%   where R_z is the Toeplitz matrix (3.15) with R(i,j) = RHO_Z^|i-j|.
%   With P = I (first choice of Section 3.3.1) this reduces to
%   A = SIGMA2 * Sigma_y^{-1}, B = Sigma_x, C = G * Sigma_x, as in
%   eq. (2.15). After solving, the mean MU removed from the entries of
%   G must be added back.
%
%   OPTS (struct with fields):
%     rho_x    - horizontal adjacent correlation (default 0.4310)
%     rho_y    - vertical adjacent correlation   (default 0.9331)
%     sigma2   - white Gaussian noise variance (default 1.825e-3)
%     psf_type - 'identity' | 'average' | 'gaussian' (default 'identity')
%     sigma_f2 - variance of the image entries; if empty it is estimated
%                from the entries of G (after removing the mean).

args = struct('rho_x', 0.4310, 'rho_y', 0.9331, ...
              'sigma2', 1.825e-3, 'psf_type', 'identity', 'sigma_f2', []);
opts = structSetDefaults(opts, args);

[N, M] = size(G);
mu = mean(G(:));
F = G - mu;

if isempty(opts.sigma_f2)
    sigma_f2 = mean(F(:).^2);
else
    sigma_f2 = opts.sigma_f2;
end

% Toeplitz matrices (3.15) and PSF (choices of Section 3.3.1)
Rx = covariance_matrix(opts.rho_x, M);
Ry = covariance_matrix(opts.rho_y, N);
Px = psf_matrix(M, opts.psf_type);
Py = psf_matrix(N, opts.psf_type);

% Sigma_z = sigma_f2 * R_z (eq. (2.8): Sigma_f = Sigma_y (kron) Sigma_x)
% The equation is written for the centered image F = G - mu (Section 3.3.1).
A = (opts.sigma2 / sigma_f2) * ((Py * Ry) \ eye(N));
B = sigma_f2 * (Px * Rx);
C = (Py \ F) * B;
end

function s = structSetDefaults(s, defaults)
%STRUCTSETDEFAULTS Adds missing fields to s with default values.
if ~isstruct(s), s = struct(); end
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i}))
        s.(f{i}) = defaults.(f{i});
    end
end
end
