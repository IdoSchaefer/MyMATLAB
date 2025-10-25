function [allfieldt, test_tpoint, propagation_grid, propagation_grid_ex_tp, i_ex_tp] = fieldw2allfieldt(fieldw, T, Nt, Nt_ts)
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    dt = T/Nt;
    t_ts = 0.5*(tcheb+1)*dt;
    tmidi = round(Nt_ts/2);
    test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
    % The time-step internal structure including the test-point:
    t_ts_tp = [t_ts(1:(Nt_ts - 1)), test_tpoint];
    allt_tp_lasti = Nt*Nt_ts + 1;
    alli = 1:allt_tp_lasti;
    i_ex_tp = alli;
    i_ex_tp(Nt_ts:Nt_ts:allt_tp_lasti) = [];
    % The propagation grid including the test point:
    propagation_grid = [kron((0:dt:(T - dt)), ones(1, Nt_ts)) + kron(ones(1, Nt), t_ts_tp), T];
    % The propagation grid excluding the test point:
    propagation_grid_ex_tp = propagation_grid(i_ex_tp);
    allfieldt = inv_dctIMtp(fieldw.', propagation_grid, T, Nt).'*sqrt(Nt*pi)/T;
end