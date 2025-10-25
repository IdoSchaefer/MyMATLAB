function [D2coef, Dcoef] = D2nschebbcoef(Nx, xdlength, alpha, beta)
% The function returns the coefficients required for the computation of the
% second derivative in a nonsymmetric modified Chebyshev grid.
    xcheb = cos((0:Nx).'*pi/Nx);
    a = 0.5*(gtild(1, alpha, beta) - gtild(-1, alpha, beta));
    D2coef = (1 - alpha*xcheb).*(1 + beta*xcheb)*a^2/(alpha*beta);
    Dcoef = (-2*alpha*beta*xcheb - alpha + beta)*a^2/(alpha*beta*xdlength);
end