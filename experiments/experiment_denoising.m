%EXPERIMENT_DENOISING Reconstruction of a noisy image via the Sylvester
% equation (Section 3.3.3 of the thesis).
%
% Compares:
%   1) the degraded image G = F + white Gaussian noise;
%   2) the full reconstruction A*Fhat + Fhat*B = C (Bartels-Stewart
%      algorithm, lyap function);
%   3) the reconstruction with truncated SVD A_k1*Fhat + Fhat*B_k2 = C
%      (Section 3.2).
%
% The parameters are those of the experiments of the thesis: 350 x 350
% image, rho_x = 0.4310, rho_y = 0.9331, noise variance
% sigma2 = 1.825e-3, P = I, and in A_k1 the last 103 singular values
% out of 350 are truncated, while nothing is truncated in B_k2.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'src'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

rng(2023);

%% Image and noise
n = 350;
img = load_example_image(n, n);
par = struct('rho_x', 0.4310, 'rho_y', 0.9331, ...
             'sigma2', 1.825e-3, 'psf_type', 'identity');

G = img + sqrt(par.sigma2) * randn(size(img));

%% Construction of A, B, C (Section 3.3.1)
[A, B, C, mu, sigma_f2] = build_sylvester_data(G, par);

%% Full solution: A*Fhat + Fhat*B = C (Bartels-Stewart)
Fhat = sylvester_solve(A, B, C);
F_full = real(Fhat) + mu;

%% Truncated SVD solution: A_k1*Fhat + Fhat*B_k2 = C
% From the thesis: 103 singular values truncated in A out of 350, none in B.
k1 = size(A, 1) - 103;
k2 = size(B, 1);
Fhat_tsvd = sylvester_tsvd_solve(A, B, C, k1, k2);
F_tsvd = Fhat_tsvd + mu;

%% Relative errors in the Frobenius norm
rel_noisy = norm(G   - img, 'fro') / norm(img, 'fro');
rel_full  = norm(F_full - img, 'fro') / norm(img, 'fro');
rel_tsvd  = norm(F_tsvd - img, 'fro') / norm(img, 'fro');
fprintf('Image variance (after mean removal): %.4e\n', sigma_f2);
fprintf('Relative error, degraded image    : %.4e\n', rel_noisy);
fprintf('Relative error, full solution    : %.4e\n', rel_full);
fprintf('Relative error, TSVD solution    : %.4e  (k1 = %d, k2 = %d)\n', ...
        rel_tsvd, k1, k2);

%% Figures (3.2, 3.3, 3.4 of the thesis)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);

subplot(2,2,1); show_image(img);    title('Original image');
subplot(2,2,2); show_image(G);      title('Image with white Gaussian noise');
subplot(2,2,3); show_image(F_full); title('Reconstructed image');
subplot(2,2,4); show_image(F_tsvd); title('Reconstructed image, truncated SVD');

saveas(fig, fullfile(outdir, 'ricostruzione.png'));

%% Figure 3.5: singular values of A and B
[U1, S1] = svd(A, 'econ');
[U2, S2] = svd(B, 'econ');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 450]);
semilogy(diag(S1), 'b', 'LineWidth', 1.2); hold on;
semilogy(diag(S2), 'r', 'LineWidth', 1.2);
grid on;
xline(k1, 'k--', sprintf('k_1 = %d', k1));
legend('\sigma_i(A)', '\sigma_i(B)', 'Location', 'southwest');
xlabel('i'); ylabel('singular values');
title('Singular values of A and B');
saveas(fig, fullfile(outdir, 'valori_singolari_AB.png'));

fprintf('Figures saved in %s\n', outdir);
close all;

function show_image(X)
%SHOW_IMAGE Displays a grayscale image.
imagesc(X); axis image off; colormap(gray);
caxis([0 1]);
end
