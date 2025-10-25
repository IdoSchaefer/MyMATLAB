function cD = Dchebc(c, domainL)
% The function computes the Chebyshev coefficients cD of the derivative of
% a signal from the Chebyshev coefficients of the signal c.
% domainL: the length of the original signal, before the mapping to the
% domain [-1, 1].
    szc = size(c);
    [N, dN] = max(szc);
    if dN == 1
        cD = zeros(N, 1);
    else
        cD = zeros(1, N);
    end
    absc = abs(c);
    c((absc/max(absc))<10*eps) = 0;
    S = zeros(1, N + 1);
    for Si = (N - 1):-1:1
        S(Si) = S(Si + 2) + Si*(c(Si + 1));
    end
    cD(1) = S(1)*2/domainL;
    cD(2:N) = S(2:N)*4/domainL;
    %cD = [S(1), 2*S(2:N)]*2/domainL;
end