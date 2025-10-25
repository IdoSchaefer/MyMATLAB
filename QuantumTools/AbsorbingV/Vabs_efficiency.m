function result = Vabs_efficiency(V, xdomain, kdomain, Nk)
% The function returns a parameter which indicates the efficiency of an
% absorbing potential. The parmeter chosen is the root mean square of the
% total norms smpled in different k values.
% V: A vector that contains the potential values at equally spaced points
% in the x domain (including the boudaries).
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nk: The number of equally spaced points in the k domain, (including the
% boundaries).
    [allT, allR, allnorm] = allTRcoef(V, xdomain, kdomain, Nk);
    result = sqrt(sum(allnorm(isfinite(allnorm)).^2));
end