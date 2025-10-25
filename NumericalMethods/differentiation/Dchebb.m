function Dv = Dchebb(v, domainL)
% The function returns the derivative of a function sampled at the
% Chebyshev points that include the boundaries of the domain.
% v: A vector that contains the function values at the Chebyshev points.
% domainL: The length of the x domain.
    N = length(v) - 1;
    cv = chebcbv(v);
    cDv = Dchebc(cv, domainL);
    cDv([1, N+1]) = 2*cDv([1, N+1]);
    Dv = dctI(sqrt(N/2)*cDv);
end