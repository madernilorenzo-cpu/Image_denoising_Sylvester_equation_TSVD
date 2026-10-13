function test_build_sylvester_data()
%TEST_BUILD_SYLVESTER_DATA Checks the construction of A, B, C.
%   1) with P = I: A = sigma2 * Sigma_y^{-1}, B = Sigma_x, C = G * Sigma_x
%      (eq. (2.15)), with Sigma_z = sigma_f2 * R_z;
%   2) for sigma2 -> 0 the solution tends to G (the filter does not
%      correct);
%   3) the matrices R_z are Toeplitz, symmetric, with unit diagonal.

rng(2);
n = 20; m = 15;
G = randn(n, m);
par = struct('rho_x', 0.4310, 'rho_y', 0.9331, ...
             'sigma2', 1e-3, 'psf_type', 'identity', 'sigma_f2', 0.5);

[A, B, C, mu] = build_sylvester_data(G, par);

Ry = covariance_matrix(par.rho_y, n);
Rx = covariance_matrix(par.rho_x, m);
assert(norm(A - (par.sigma2/par.sigma_f2)*inv(Ry), 'fro')/norm(A,'fro') < 1e-12);
assert(norm(B - par.sigma_f2*Rx, 'fro')/norm(B,'fro') < 1e-12);
assert(norm(C - (G - mean(G(:)))*(par.sigma_f2*Rx), 'fro')/norm(C,'fro') < 1e-10);
assert(abs(mu - mean(G(:))) < 1e-15);

% R_z: Toeplitz, unit diagonal, symmetric
for rho = [0.4310 0.9331]
    R = covariance_matrix(rho, 8);
    assert(all(diag(R) == 1));
    assert(isequal(R, R.'));
    assert(isequal(R, toeplitz(R(:,1))));
end

% sigma2 -> 0: the solution goes back to G
As = A * 1e-12;
Fh = sylvester_solve(As, B, C);
assert(norm(Fh + mu - G, 'fro') / norm(G, 'fro') < 1e-8);

fprintf('test_build_sylvester_data: OK\n');
end
