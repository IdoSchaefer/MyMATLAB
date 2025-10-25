function w = modchebbintegw(Nx, p, xdlength)
% The function computes the integration weights w for a modified Chebyshev
% grid.
% Nx: The number of grid points minus 1.
% p: The streching parameter from the new article, 0<p<pi/2
% xdlength: The x grid length
    theta  = (0:pi/Nx:pi).';
    w = sin(p)/p*xdlength/2*sin(theta)./sqrt(1 - cos(theta).^2*sin(p)^2)*pi/Nx;
end