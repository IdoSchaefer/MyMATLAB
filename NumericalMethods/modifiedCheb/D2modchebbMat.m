function D2M = D2modchebbMat(Nx, xdlength, p)
% The function computes the second derivation matrix D2M for a boundary
% including Chebyshev grid.
% Nx: The number of grid points minus 1
% xdlength: The length of the x grid
% p: The streching parameter, 0<p<pi/2
    xcheb = cos((0:Nx).'*pi/Nx);
    % The inverse of dx/d(xcheb): 
    Dxinv = p*sqrt(1-(sin(p)*xcheb).^2)/sin(p); 
    D2M = diag(Dxinv.^2)*D2chebbMat(Nx, xdlength) - ...
        p^2*diag(xcheb)*DchebbMat(Nx, xdlength)*2/xdlength;
end