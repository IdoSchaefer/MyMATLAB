function [DM, DMcheb] = DchebbMat(Nx, xdlength)
% The function computes the first derivation matrix in a Chebyshev grid
% that includes the boundaries.
% Input:
% Nx: number of grid points minus 1.
% xdlength: the length of the x grid
% Output:
% DM: The derivation matrix in the x domain
% DMcheb: the derivation matrix in the Chebyshev domain
    l = zeros(1, Nx + 1);
    l(2:(Nx + 1)) = (1:Nx)*2;
    DMcheb = zeros(Nx + 1, Nx + 1);
    for ni = 1:(Nx + 1)
        DMcheb(ni, (ni + 1):2:(Nx + 1)) = l((ni + 1):2:(Nx + 1));
    end
    DMcheb(1, :) = DMcheb(1, :)/2;
    DMcheb = DMcheb*2/xdlength;
    DM = dctIM([2*DMcheb(1, :);
                DMcheb(2:Nx, :);
                2*DMcheb(Nx + 1, :)])';
    DM = dctIM([DM(1, :);
                DM(2:Nx, :);
                DM(Nx + 1, :)])';
    DM(:, [1, Nx + 1]) = DM(:, [1, Nx + 1])/2;
end