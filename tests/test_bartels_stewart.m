function test_bartels_stewart()
%TEST_BARTELS_STEWART Checks the Bartels-Stewart algorithm.
%   1) against the vectorized solution with H = kron(I,A)+kron(B',I);
%   2) against lyap (if the Control System Toolbox is installed);
%   3) residual of the equation A*X + X*B = C.

rng(0);
n = 30; m = 25;

% A and B symmetric positive definite: lambda_i(A) + lambda_j(B) > 0,
% so the equation has a unique solution.
A = randn(n); A = (A + A')/2 + n*eye(n);
B = randn(m); B = (B + B')/2 + m*eye(m);
C = randn(n, m);

X  = bartels_stewart(A, B, C);

% Reference: vectorized solution (H\c) with H = kron(I,A)+kron(B',I)
H = kron_operator(A, B);
Xv = reshape(H \ C(:), n, m);
err_kron = norm(X - Xv, 'fro') / norm(Xv, 'fro');
assert(err_kron < 1e-10, 'Mismatch with H\\c: %g', err_kron);

% If lyap (Control System Toolbox) is available, compare against it too.
if exist('lyap', 'file')
    Xr = lyap(A, B, -C);
    err_lyap = norm(X - Xr, 'fro') / norm(Xr, 'fro');
    assert(err_lyap < 1e-10, 'Mismatch with lyap: %g', err_lyap);
else
    err_lyap = NaN;
end

res = norm(A*X + X*B - C, 'fro') / norm(C, 'fro');
assert(res < 1e-10, 'Residual too large: %g', res);

fprintf('test_bartels_stewart: OK (err lyap %.2e, err kron %.2e, res %.2e)\n', ...
        err_lyap, err_kron, res);
end
