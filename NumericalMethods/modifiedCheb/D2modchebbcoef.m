function [D2coef, Dcoef] = D2modchebbcoef(Nx, xdlength, p)
% The function returns the coefficients required for the computation of the
% second derivative in a modified Chebyshev grid.
    xcheb = cos((0:Nx).'*pi/Nx);
%     % The inverse of dx/d(xcheb): 
%     Dxinv = p*sqrt(1-(sin(p)*xcheb).^2)/sin(p);
%     D2coef = Dxinv.^2;
% This is better:
    D2coef = (p/sin(p))^2*(1-(sin(p)*xcheb).^2);
    Dcoef = -p^2*xcheb*2/xdlength;
end