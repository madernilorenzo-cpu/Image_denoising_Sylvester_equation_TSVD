function X = bartels_stewart(A, B, C)
%BARTELS_STEWART Solves the Sylvester equation A*X + X*B = C.
%
%   X = BARTELS_STEWART(A, B, C) implements the Bartels-Stewart
%   algorithm (Section 2.1): A and B are reduced to Schur form, the
%   transformed equation T1*Y + Y*T2 = C~ is solved by substitution and
%   the result is finally transformed back.
%
%   Here the complex Schur form is used (always upper triangular); the
%   original algorithm uses the real Schur form with 2x2 blocks and
%   linear systems of order at most four, see eq. (2.7). MATLAB's lyap
%   implements the real variant.
%
%   The equation has a unique solution if lambda_i(A) + lambda_j(B) ~= 0
%   for every pair (i, j).

[n, n2] = size(A);
[m, m2] = size(B);
if n ~= n2 || m ~= m2 || size(C, 1) ~= n || size(C, 2) ~= m
    error('bartels_stewart:wrongDims', ...
        'Need A (n x n), B (m x m), C (n x m).');
end

% Complex Schur form: A = U*T1*U', B = V*T2*V' with T1, T2 upper
% triangular.
[U, T1] = schur(A, 'complex');
[V, T2] = schur(B, 'complex');

% C~ = U' * C * V
Ct = U' * C * V;

% Column-wise substitution (eq. (2.7) in the triangular case): the j-th
% column of Y satisfies
%   (T1 + T2(j,j)*I) * Y(:,j) = C~(:,j) - Y(:,1:j-1) * T2(1:j-1,j)
Y = zeros(n, m, 'like', Ct);
In = eye(n);
for j = 1:m
    rhs = Ct(:, j);
    if j > 1
        rhs = rhs - Y(:, 1:j-1) * T2(1:j-1, j);
    end
    Mj = T1 + T2(j, j) * In;
    if min(abs(diag(T1) + T2(j, j))) < eps * norm(T1, 'fro')
        warning('bartels_stewart:nearlySingular', ...
            'T1 + T2(%d,%d)*I is nearly singular.', j, j);
    end
    Y(:, j) = Mj \ rhs;
end

% Back-transform: X = U * Y * V'
X = U * Y * V';
end
