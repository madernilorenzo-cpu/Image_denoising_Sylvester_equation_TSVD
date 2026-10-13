%EXPERIMENT_APPROXH Approximation of the operator H = I (kron) A + B' (kron) I
% (Section 3.3.2 of the thesis).
%
% Compares the decay of the singular values of the operator H with the
% vector ThetaA+B of all pairwise sums of singular values of A and B
% (sorted in decreasing order):
%
%   SigmaH    = svd(H),        H = kron(I, A) + kron(B', I)
%   ThetaA+B  = sort( sigma_i(A) + sigma_j(B), 'descend' )
%
% By inequality (3.16) we should have
%
%   SigmaH(k) <= ThetaA+B(k)    for every k.
%
% Note: in the thesis the experiment is run at full size (350 x 350,
% i.e. H of order 122500); here a smaller image is used so that the SVD
% of H is affordable on a laptop.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'src'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

rng(7);

%% Reduced-size test problem
n = 48;   % rows    -> A is n x n
m = 40;   % columns -> B is m x m
img = load_example_image(n, m);
par = struct('rho_x', 0.4310, 'rho_y', 0.9331, ...
             'sigma2', 1.825e-3, 'psf_type', 'identity');
G = img + sqrt(par.sigma2) * randn(size(img));

[A, B, ~] = build_sylvester_data(G, par);

%% Operator H and its spectrum
H = kron_operator(A, B);
SigmaH = svd(H);

%% Sums of singular values of A and B (ThetaA+B)
sa = svd(A);
sb = svd(B);
ThetaAplusB = sort(reshape(sa(:) + sb(:).', [], 1), 'descend');

%% Check of inequality (3.16) and relative difference
assert(all(SigmaH <= ThetaAplusB + 1e-12 * max(SigmaH)), ...
    'Inequality SigmaH(k) <= ThetaA+B(k) violated.');
rel = norm(SigmaH - ThetaAplusB) / max(SigmaH);
fprintf('norm(SigmaH - ThetaA+B)/max(SigmaH) = %.4e\n', rel);

%% Plot (Figure 3.1 of the thesis)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 700 450]);
semilogy(SigmaH, 'b', 'LineWidth', 1.2); hold on;
semilogy(ThetaAplusB, 'r--', 'LineWidth', 1.2);
grid on; legend('\SigmaH', '\Theta_{A+B}', 'Location', 'southwest');
xlabel('k'); ylabel('singular values');
title('Singular values of H and sums of singular values of A and B');
saveas(fig, fullfile(outdir, 'sigmaH_vs_thetaAplusB.png'));

fprintf('Size of H: %d x %d\n', size(H, 1), size(H, 2));
close all;
