function X = sylvester_solve(A, B, C)
%SYLVESTER_SOLVE Solves the Sylvester equation A*X + X*B = C.
%
%   X = SYLVESTER_SOLVE(A, B, C) uses MATLAB's lyap function (Control
%   System Toolbox) if available, as in the experiments of the thesis;
%   otherwise it falls back to the standalone implementation of the
%   Bartels-Stewart algorithm (bartels_stewart.m), which is the
%   algorithm internally used by lyap.

if exist('lyap', 'file')
    X = lyap(A, B, -C);
else
    X = bartels_stewart(A, B, C);
end
end
