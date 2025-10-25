function D2M = D2nschebbMat(Nx, xdlength, alpha, beta)
% The function computes the second derivation matrix D2M for a boundary
% including nonsymmetric modified Chebyshev grid.
% Nx: The number of grid points minus 1
% xdlength: The length of the x grid
% alpha: The streching parameter for the right boundary
% beta: The streching parameter for the left boundary
    [D2coef, Dcoef] = D2nschebbcoef(Nx, xdlength, alpha, beta);
    D2M = diag(D2coef)*D2chebbMat(Nx, xdlength) + ...
        diag(Dcoef)*DchebbMat(Nx, xdlength);
end