function D2v = D2modchebb(v, domainL, p)
% The function returns the second derivative of a function sampled at the
% modified Chebyshev points that include the bounderies.
% v: A vector that contains the function values at the Chebyshev points.
% domainL: The length of the x domain.
    Nx = length(v) - 1;
    xcheb = cos((0:Nx).'*pi/Nx);
    % The inverse of dx/d(xcheb): 
    Dxinv = p*sqrt(1-(sin(p)*xcheb).^2)/sin(p); 
    D2v = Dxinv.^2.*D2chebb(v, domainL) - ...
        p^2*xcheb.*Dchebb(v, domainL)*2/domainL;    
end