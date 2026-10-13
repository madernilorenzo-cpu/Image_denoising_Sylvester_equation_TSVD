function test_sylvester_tsvd()
%TEST_SYLVESTER_TSVD Checks the truncated-SVD solution (Section 3.2).
%   1) without truncation (k1, k2 full) it matches the exact solution;
%   2) truncating only A, the equation A_k1*F + F*B_k2 = C stays
%      solvable and the residual is zero;
%   3) as k1 grows the solution converges to the one of the full
%      equation;
%   4) truncating both A and B makes the problem singular (eigenvalues
%      0+0): the function must report it (Section 3.3.3 of the thesis).

rng(1);
n = 40; m = 32;
A = randn(n); A = (A + A')/2 + n*eye(n);
B = randn(m); B = (B + B')/2 + m*eye(m);
C = randn(n, m);

% 1) no truncation -> exact solution
Xfull = sylvester_tsvd_solve(A, B, C, n, m);
Xr = sylvester_solve(A, B, C);
err = norm(Xfull - Xr, 'fro') / norm(Xr, 'fro');
assert(err < 1e-10, 'With full k it must match the exact solution: %g', err);

% 2) truncation of A only: zero residual
[U1, S1, V1] = svd(A, 'econ');
prev_gap = inf;
for k1 = [n/2, 3*n/4, n-4]
    Xk = sylvester_tsvd_solve(A, B, C, k1, m);
    Ak1 = U1(:,1:k1) * S1(1:k1,1:k1) * V1(:,1:k1)';
    res = norm(Ak1*Xk + Xk*B - C, 'fro') / norm(C, 'fro');
    assert(res < 1e-10, 'Residual of the truncated equation not zero: %g', res);

    % 3) convergence to the full solution as k1 grows
    gap = norm(Xk - Xr, 'fro') / norm(Xr, 'fro');
    fprintf('  k1 = %2d: residual %.2e, distance from exact solution %.3e\n', ...
            k1, res, gap);
    assert(gap < prev_gap, 'Distance from the solution does not decrease with k1.');
    prev_gap = gap;
end

% 4) simultaneous truncation of A and B: singular problem
try
    sylvester_tsvd_solve(A, B, C, n-1, m-1);
    assert(false, 'Expected error for simultaneous truncation of A and B.');
catch ME
    assert(contains(ME.identifier, 'singularProblem'), ...
        'Unexpected error: %s', ME.message);
end

fprintf('test_sylvester_tsvd: OK\n');
end
