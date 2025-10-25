function [chis, Jmax] = chis_gate(psis, targets, limits)
    N = length(psis);
    chis = zeros(N, 1);
    lower_limit = 1;
    overlaps = zeros(1, 4);
    for ipsi = 1:3
        upper_limit = limits(ipsi);
        overlaps(ipsi) = targets(lower_limit:upper_limit)'*psis(lower_limit:upper_limit);
        lower_limit = upper_limit + 1;
    end
    upper_limit = N;
    overlaps(4) = targets(lower_limit:upper_limit)'*psis(lower_limit:upper_limit);
    chis(1:limits(1)) = overlaps(2)*overlaps(3)*conj(overlaps(4))*targets(1:limits(1))/4;
    chis((limits(1) + 1):limits(2)) = overlaps(1)*overlaps(4)*conj(overlaps(3))*targets((limits(1) + 1):limits(2))/4;    
    chis((limits(2) + 1):limits(3)) = overlaps(1)*overlaps(4)*conj(overlaps(2))*targets((limits(2) + 1):limits(3))/4;    
    chis((limits(3) + 1):N) = overlaps(2)*overlaps(3)*conj(overlaps(1))*targets((limits(3) + 1):N)/4;
    Jmax = 0.5 + real(overlaps(1)*overlaps(4)*conj(overlaps(2)*overlaps(3)))/2;
end