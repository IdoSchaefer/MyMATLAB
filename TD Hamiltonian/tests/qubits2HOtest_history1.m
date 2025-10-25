load qubits2HOtest_problem
Nt_ts = 9;
tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
dt = 1e-2*2*pi;
t_ts = 0.5*(tcheb+1)*dt;
tmidi = round(Nt_ts/2);
test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
% The time-step internal structure including the test-point:
t_ts_tp = [t_ts(1:(Nt_ts - 1)), test_tpoint];
Nt = T2mb2/dt;
dctfactor = T2mb2/(sqrt(Nt*pi));
allt_tp_lasti = Nt*Nt_ts + 1;
allfield = zeros(2, allt_tp_lasti);
for icontrol = 1:2
    allfield(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t_ts_tp)/dctfactor;
end
propagation_grid = [kron((0:dt:(T2mb2 - dt)), ones(1, Nt_ts)) + kron(ones(1, Nt), [t_ts(1:(Nt_ts - 1)), test_tpoint]), T2mb2];
i_without_tp = 1:allt_tp_lasti;
i_without_tp(Nt_ts:Nt_ts:allt_tp_lasti) = [];
[U, mniterc, matvecs, max_errors, history] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 2/3, allfield, [], [-25 35], u02modes, [0, T2mb2], Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
t01tp = [t01(1:(Nt_ts-1)), (t01(tmidi) + t01(tmidi + 1))/2];
allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(140*9 + 1))*2, T2mb2, 140).'*sqrt(140*pi)/T2mb2;
propagation_grid_max = propagation_grid(1:(140*9 + 1))*2;
plot(propagation_grid(i_without_tp)/(2*pi), allfield(:,i_without_tp))
hold on
plot(propagation_grid_max(i_without_tp(1:1121))/(2*pi), allfieldmax(:,i_without_tp(1:1121)))
[Umax, mnitercmax, matvecsmax, max_errorsmax, historymax] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfieldmax, [], [-15 35], u02modes, [0, T2mb2], Nt/2, 9, 9, 1e-6, 10, 16, test_tpoint*2, false);
allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(70*9 + 1))*4, T2mb2, 70).'*sqrt(70*pi)/T2mb2;