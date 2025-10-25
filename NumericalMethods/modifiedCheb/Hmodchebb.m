function H = Hmodchebb(Vf, domainL, Nx, p, m)
% The function computes the Hamoltonian matrix for the modified Chebyshev points.
% Vf: The potential energy function. It's a function handle.
% domainL: The length of the x domain.
% m: The mass
    x = getmodchebbp(Nx, p, -domainL/2, domainL/2);
    H = diag(Vf(x)) - D2modchebbMat(Nx, domainL, p)/(2*m);
end