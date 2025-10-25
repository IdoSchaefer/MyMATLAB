function [parameter, gradient] = Vabs_efficiency1(V, xdomain, kdomain, Nk)
% The function returns a parameter which indicates the efficiency of an
% absorbing potential. The parmeter chosen is the sum of the
% squares of the norms sampled in different k values.
% V: A vector that contains the potential values at equally spaced points
% in the x domain (including the boudaries).
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nk: The number of equally spaced points in the k domain, (including the
% boundaries).
    if nargout>1
        [allT, allR, allnorm, allgradnorm] = allTRcoef1(V, xdomain, kdomain, Nk);
        gradient = 2*sum(allgradnorm*diag(allnorm), 2);
    else
        [allT, allR, allnorm] = allTRcoef1(V, xdomain, kdomain, Nk);
    end    
    % parameter = sum(allnorm(isfinite(allnorm)).^2);
    parameter = sum(allnorm.^2);
end