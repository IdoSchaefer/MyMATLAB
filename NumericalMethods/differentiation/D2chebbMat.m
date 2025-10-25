function [D2M, D2Mcheb] = D2chebbMat(Nx, xdlength)
% The function computes the second derivation matrix in a Chebyshev grid
% that includes the boundaries.
% Input:
% Nx: number of grid points
% xdlength: the length of the x grid
% Output:
% D2M: The second derivation matrix in the x domain
% D2Mcheb: the second derivation matrix in the Chebyshev domain
    l = zeros(1, Nx + 1);
    l(3:(Nx + 1)) = 2:Nx;
    D2Mcheb = zeros(Nx + 1, Nx + 1);
    for ni = 1:(Nx + 1)
        D2Mcheb(ni, (ni + 2):2:(Nx + 1)) = l((ni + 2):2:(Nx + 1)).*(l((ni + 2):2:(Nx + 1)).^2 - (ni - 1)^2);
    end
    D2Mcheb(1, :) = D2Mcheb(1, :)/2;
    D2Mcheb = D2Mcheb*4/xdlength^2;
    D2M = dctIM([2*D2Mcheb(1, :);
                 D2Mcheb(2:Nx, :);
                 2*D2Mcheb(Nx + 1, :)])';
    D2M = dctIM([D2M(1, :);
                 D2M(2:Nx, :);
                 D2M(Nx + 1, :)])';
    D2M(:, [1, Nx + 1]) = D2M(:, [1, Nx + 1])/2;
end