function D2v = D2chebb(v, domainL)
% The function returns the second derivative of a function sampled at the
% Chebyshev points that include the bounderies.
% v: A vector that contains the function values at the Chebyshev points.
% domainL: The length of the x domain.
    N = length(v) - 1;
    cv = chebcbv(v);
    cD2v = D2chebc(cv, domainL);
    cD2v([1, N+1]) = 2*cD2v([1, N+1]);
    D2v = dctI(sqrt(N/2)*cD2v);
end