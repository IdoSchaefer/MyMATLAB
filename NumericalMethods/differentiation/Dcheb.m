function Dv = Dcheb(v, domainL)
% The function returns the derivative of a function sampled at the
% Chebyshev points.
% v: A vector that contains the function values at the Chebyshev points.
% domainL: The length of the x domain.
    N = length(v);
    cv = chebcv(v);
    cDv = Dchebc(cv, domainL);
    cDv(1) = cDv(1)*sqrt(N);
    cDv(2:N) = cDv(2:N)*sqrt(N/2);
    Dv = idct(cDv);
    %Dv = idct([cDv(1), cDv(2:N)/sqrt(2)]*sqrt(N));
end