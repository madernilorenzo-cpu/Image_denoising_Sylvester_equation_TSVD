function H = kron_operator(A, B)
%KRON_OPERATOR Operator H = I (kron) A + B' (kron) I of the thesis.
%
%   H = KRON_OPERATOR(A, B) builds the matrix
%
%       H = I_m (kron) A + B' (kron) I_n
%
%   such that, by Proposition 1.2.2 (vec(A X B) = (B' (kron) A) vec(X)),
%
%       vec(A*F + F*B) = H * vec(F),
%
%   where A is n x n, B is m x m and vec is the column-wise vectorization.
%   This is the operator of the problem g = H*f of Chapter 3: large
%   (n*m x n*m) and ill-conditioned.

[n, n2] = size(A);
[m, m2] = size(B);
if n ~= n2 || m ~= m2
    error('kron_operator:wrongDims', 'A and B must be square.');
end
H = kron(eye(m), A) + kron(B.', eye(n));
end
