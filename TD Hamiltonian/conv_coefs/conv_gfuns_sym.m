function gfuns = conv_gfuns_sym(M, p_max, xivals)
% The function computes the values of g^{(p)}(xi) at values xivals, for all
% p up to the maximal iteration number p_max. The computation is based on a
% symbolic integration. It was found that in high orders this computation
% is more stable than the iterative scheme of conv_gfuns.m. However, this
% is much more time-consuming, and the typically required orders are much less than
% those with stability issues.
% Input:
% M: The M value
% p_max: The maximal iteration number
% xivals: A row vector of the xi values to be computed
    gfuns = zeros(p_max + 1, length(xivals));
    cheb_xi = (-cos(0:pi/(M - 1):pi) + 1)/2;
    syms xi
    RMelements = xi - cheb_xi;
    RM = prod(RMelements);
    gp_sym = int(RM, xi, 0, 1 + xi);
    gfuns(1, :) = double(subs(gp_sym, xivals));
    for p_i = 2:(p_max + 1)
        gp_sym = int(gp_sym*(xi - 1/2), xi, 0, xi);
        gfuns(p_i, :) = double(subs(gp_sym, xivals));
    end
end