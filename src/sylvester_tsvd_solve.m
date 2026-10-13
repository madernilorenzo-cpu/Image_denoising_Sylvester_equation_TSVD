function X = sylvester_tsvd_solve(A, B, C, k1, k2)
%SYLVESTER_TSVD_SOLVE Solves A_k1*X + X*B_k2 = C via TSVD (Section 3.2).
%
%   X = SYLVESTER_TSVD_SOLVE(A, B, C, K1, K2) truncates the SVDs of A and
%   B to the first K1 and K2 singular values,
%
%       A_k1 = U1(:,1:K1) * S1(1:K1,1:K1) * V1(:,1:K1)'
%       B_k2 = U2(:,1:K2) * S2(1:K2,1:K2) * V2(:,1:K2)'
%
%   and solves the approximate Sylvester equation
%
%       A_k1 * X + X * B_k2 = C
%
%   with the Bartels-Stewart algorithm (lyap if available), as in the
%   experiments of the thesis.
%
%   NOTE (Section 3.3.3 of the thesis): the equation has a unique
%   solution only if lambda_i(A_k1) + lambda_j(B_k2) ~= 0 for every
%   pair. If BOTH matrices are truncated, the null eigenvalues of A_k1
%   and B_k2 sum to zero and the problem becomes singular: this is why
%   in the experiments of the thesis only A is truncated (103 out of
%   350 singular values) and no truncation of B is possible.

[U1, S1, V1] = svd(A, 'econ');
[U2, S2, V2] = svd(B, 'econ');

n = size(A, 1);
m = size(B, 1);
if k1 > n || k2 > m
    error('sylvester_tsvd_solve:kTooLarge', ...
        'k1 and k2 cannot exceed the sizes of A and B.');
end
if k1 < n && k2 < m
    error('sylvester_tsvd_solve:singularProblem', ...
        ['Truncating both A and B gives lambda_i(A_k1) + lambda_j(B_k2) = 0 ', ...
         'for the null eigenvalues: the problem is singular. ', ...
         'Truncate only one of the two matrices (see Section 3.3.3 of the thesis).']);
end

Ak1 = U1(:, 1:k1) * S1(1:k1, 1:k1) * V1(:, 1:k1)';
Bk2 = U2(:, 1:k2) * S2(1:k2, 1:k2) * V2(:, 1:k2)';

X = sylvester_solve(Ak1, Bk2, C);

if isreal(A) && isreal(B) && isreal(C)
    X = real(X);
end
end
