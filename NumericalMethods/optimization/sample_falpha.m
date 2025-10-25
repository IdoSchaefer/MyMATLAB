function [falpha, Dfalpha] = sample_falpha(fDf, x0, direction, alpha)
    Nalpha = length(alpha);
    falpha = zeros(1, Nalpha);
    Dfalpha = zeros(1, Nalpha);
    for alphai = 1:Nalpha
        xalpha = x0 + alpha(alphai)*direction;
        [falpha(alphai), fgrad_xalpha] = fDf(xalpha);
        Dfalpha(alphai) = fgrad_xalpha.'*direction;
    end
    figure
    plot(alpha, falpha)
    figure
    plot(alpha, Dfalpha)
end