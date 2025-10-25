function cD2 = D2chebc(c, domainL)
% The function computes the Chebyshev coefficients cD2 of the second derivative of
% a signal from the Chebyshev coefficients of the signal c.
% domainL: the length of the original signal, before the mapping to the
% domain [-1, 1].
    szc = size(c);
    [N, dN] = max(szc);
    if dN == 1
        cD2 = zeros(N + 2, 1);
    else
        cD2 = zeros(1, N + 2);
    end
    absc = abs(c);
    c((absc/max(absc))<10*eps) = 0;
    for cD2i = (N - 2):-1:1
        cD2(cD2i) = 4*(cD2i + 1)*cD2i*c(cD2i + 2) + ...
            (2*(cD2i + 1)*cD2(cD2i + 2) - cD2i*cD2(cD2i + 4))/(cD2i + 2);
    end
    cD2(1) = cD2(1)/2;
    cD2 = cD2(1:N)*4/domainL^2;
    %cD2 = [cD2(1)/2, cD2(2:N)]*4/domainL^2;
end
    