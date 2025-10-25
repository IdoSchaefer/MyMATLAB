function H = Hnschebb(Vf, domainL, Nx, alpha, beta, m)
% The function computes the Hamoltonian matrix for the modified Chebyshev points.
% Vf: The potential energy function. It's a function handle.
% domainL: The length of the x domain.
% m: The mass
    x = getnschebbp(Nx, alpha, beta, -domainL/2, domainL/2);
    H = diag(Vf(x)) - D2nschebbMat(Nx, domainL, alpha, beta)/(2*m);
end