function [D2M, D2Mcheb] = D2chebMat(Nx, xdlength)
% The function computes the second derivation matrix in a Chebyshev grid.
% Input:
% Nx: number of grid points
% xdlength: the length of the x grid
% Output:
% D2M: The second derivation matrix in the x domain
% D2Mcheb: the second derivation matrix in the Chebyshev domain
    l = zeros(1, Nx);
    l(3:Nx) = 2:(Nx - 1);
    D2Mcheb = zeros(Nx, Nx);
    for ni = 1:Nx
        D2Mcheb(ni, (ni + 2):2:Nx) = l((ni + 2):2:Nx).*(l((ni + 2):2:Nx).^2 - (ni - 1)^2);
    end
    D2Mcheb(1, :) = D2Mcheb(1, :)/2;
    D2Mcheb = D2Mcheb*4/xdlength^2;
    D2M = idct([D2Mcheb(1, :);
                D2Mcheb(2:Nx, :)/sqrt(2)]*sqrt(Nx))';
    D2M = idct([D2M(1, :);
                D2M(2:Nx, :)*sqrt(2)]/sqrt(Nx))';
end