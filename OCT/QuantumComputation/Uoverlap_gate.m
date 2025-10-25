function [Jmax_overlap, phases] = Uoverlap_gate(U, targets, limits)
    phases = zeros(4, 1);
    N = length(U);
    targets_phases = targets;
    lower_limit = 1;
    for ipsi = 1:3
        upper_limit = limits(ipsi);
        indices_psii = lower_limit:upper_limit;
        phases(ipsi) = targets(indices_psii)'*U(indices_psii);
        phases(ipsi) = phases(ipsi)/abs(phases(ipsi));
        targets_phases(indices_psii) = targets(indices_psii)*phases(ipsi);
        lower_limit = upper_limit + 1;
    end
    upper_limit = N;
    indices_psii = lower_limit:upper_limit;
    phases(4) = targets(indices_psii)'*U(indices_psii);
    phases(4) = phases(4)/abs(phases(4));
    targets_phases(indices_psii) = targets(indices_psii)*phases(4);
    Jmax_overlap = targets_phases'*U/4;
    Jmax_overlap = Jmax_overlap*conj(Jmax_overlap);
end