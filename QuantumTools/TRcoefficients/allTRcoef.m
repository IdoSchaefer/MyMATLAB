function [allT, allR, allnorm] = allTRcoef(V, xdomain, kdomain, Nk)
% The function returns the transmission and reflection coefficients for a
% domain of wave numbers, sampled at equally spaced points.
% Input:
% V: A vector that contains the potential values at equally spaced points
% in the x domain (including the boudaries).
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nk: The number of equally spaced points in the k domain, (including the
% boundaries).
% Output:
% allT: A vector that contains all the transmission coefficients in all k.
% allR: A vector that contains all the reflection coefficients in all k.
% allnorm: A vector that contains all the total norms in all k.
    Nx = length(V);
    dk = (kdomain(2) - kdomain(1))/(Nk - 1);
    k = kdomain(1):dk:kdomain(2);
    allT = zeros(1, Nk);
    allR = zeros(1, Nk);
    allnorm = zeros(1, Nk);
    for ki = 1:Nk
        [allT(ki), allR(ki), allnorm(ki)] = TRcoef(k(ki)^2/2, V, xdomain, Nx);
    end
end