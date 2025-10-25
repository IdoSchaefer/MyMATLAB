function [chiihterm, Jforb] = fJforb(allpsi, Nt_ts, integw, penalforbv)
    allt_tp_lasti = size(allpsi, 2);
    Nt = (allt_tp_lasti - 1)/Nt_ts;
    % The time indices of the propagation grid excluding the test-points:
    indices_excluding_tp = 1:allt_tp_lasti;
    indices_excluding_tp(Nt_ts:Nt_ts:Nt*Nt_ts) = [];
    chiihterm = penalforbv.*allpsi;
    exvals = real(sum(conj(allpsi(:, indices_excluding_tp)).*chiihterm(:, indices_excluding_tp)));
    Jforb = -sum(integw.*exvals);
end