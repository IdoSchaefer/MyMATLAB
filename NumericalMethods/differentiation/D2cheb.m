function D2v = D2cheb(v, domainL)
% The function returns the second derivative of a function sampled at the
% Chebyshev points.
% v: A vector that contains the function values at the Chebyshev points.
% domainL: The length of the x domain.
    N = length(v);
    cv = chebcv(v);
    cD2v = D2chebc(cv, domainL);
    cD2v(1) = cD2v(1)*sqrt(N);
    cD2v(2:N) = cD2v(2:N)*sqrt(N/2);
    D2v = idct(cD2v);
    %D2v = idct([cD2v(1), cD2v(2:N)/sqrt(2)]*sqrt(N));
end