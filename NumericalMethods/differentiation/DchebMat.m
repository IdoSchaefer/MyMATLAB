function [DM, DMcheb] = DchebMat(Nx, xdlength)
% The function computes the first derivation matrix in a Chebyshev grid.
% Input:
% Nx: number of grid points
% xdlength: the length of the x grid
% Output:
% DM: The derivation matrix in the x domain
% DMcheb: the derivation matrix in the Chebyshev domain
    l = zeros(1, Nx);
    l(2:Nx) = (1:(Nx - 1))*2;
    DMcheb = zeros(Nx, Nx);
    for ni = 1:Nx
        DMcheb(ni, (ni + 1):2:Nx) = l((ni + 1):2:Nx);
    end
    DMcheb(1, :) = DMcheb(1, :)/2;
    DMcheb = DMcheb*2/xdlength;
    DM = idct([DMcheb(1, :);
                DMcheb(2:Nx, :)/sqrt(2)]*sqrt(Nx))';
    DM = idct([DM(1, :);
                DM(2:Nx, :)*sqrt(2)]/sqrt(Nx))';
end