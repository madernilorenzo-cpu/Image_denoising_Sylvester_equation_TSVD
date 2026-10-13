function R = covariance_matrix(rho, n)
%COVARIANCE_MATRIX Toeplitz matrix R_z of the thesis, eq. (3.15).
%
%   R = COVARIANCE_MATRIX(RHO, N) returns the N x N Toeplitz matrix
%   defined in Section 3.3.1:
%
%       R(i,j) = RHO^|i-j|,   R(i,i) = 1.
%
%   RHO is the adjacent-entry correlation (rho_x or rho_y in the
%   thesis). For |RHO| < 1 the matrix is positive definite.

v = rho.^(0:n-1);
R = toeplitz(v);
end
