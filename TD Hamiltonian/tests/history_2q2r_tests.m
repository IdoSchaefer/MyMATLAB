0.25/0.03
load HTLSs_harmonic
g2GHz
g1GHz
plot(0:0.01:2.8, fieldt2mb2)
0.2/g2GHz
Delta_nuGHz/g2GHz
H0
mu
[Gop, Gdiff_op] = Hmats2Gop_dipole(H0, mu, @(u,t) cos(t))
Gop([1;0], pi, [1;0])
H0*[1;0] - cos(pi)*mu*[1;0]
Gdiff_op([1,0; 0, 1], [pi 3*pi/2], 1/sqrt(2)*[1;1], 0)
doc feval
clear Gop Gdiff_op
plot(0:0.01:2.8, fieldt2mb2)
u02modes
save qubits2HOtest_problem fieldw2mb2 u02modes H2modesa H1_2modesa H2_2modesa Hu2modesa H1u2modesa H2u2modesa Hoperations2modesa g2GHz Delta_nuGHz
T
T2mb2 = 2.8*2*pi
save qubits2HOtest_problem fieldw2mb2 u02modes H2modesa H1_2modesa H2_2modesa Hu2modesa H1u2modesa H2u2modesa Hoperations2modesa g2GHz Delta_nuGHz T2mb2
%-- 19/07/2022 9:32 --%
load qubits2HOtest_problem
whos
u02modes
options_gate4
load HTLSs_harmonic
options_gate4
clear all
load qubits2HOtest_problem
Nt_ts = 9;
tcheb = -cos(((1:Nt_ts).' - 1)*pi/(Nt_ts-1));
dt = 1e-2*2*pi;
t_ts = 0.5*(tcheb+1)*dt;
tmidi = round(Nt_ts/2);
test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
% The time-step internal structure including the test-point:
t_ts_tp = [t_ts(1:(Nt_ts - 1)); test_tpoint];
Nt = T2mb2/dt
dctfactor = T2mb2/(sqrt(Nt*pi));
allt_tp_lasti = Nt*Nt_ts + 1;
allfield = zeros(2, allt_tp_lasti);
for icontrol = 1:2, allfield(icontrol, :) = dctIintgrid1(fieldw2mb2(:, icontrol).', T2mb2, t_ts_tp.')/dctfactor; end
allt_tp_lasti
for icontrol = 1:2, allfield(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t_ts_tp.')/dctfactor; end
t_ts_tp
propagation_grid = [kron((0:dt:(T2mb2 - dt)), ones(1, Nt_ts)) + kron(ones(1, Nts), [t_ts(1:(Nt_ts - 1)).', test_tpoint]), T2mb2];
propagation_grid = [kron((0:dt:(T2mb2 - dt)), ones(1, Nt_ts)) + kron(ones(1, Nt), [t_ts(1:(Nt_ts - 1)).', test_tpoint]), T2mb2];
i_without_tp = 1:all_tp_lasti;
i_without_tp = 1:allt_tp_lasti;
i_without_tp(Nt_ts:Nt_ts:allt_tp_lasti) = [];
i_without_tp(1:30)
propagation_grid(i_without_tp(1:30))
figure
plot(propagation_grid(i_without_tp), allfield(:, i_without_tp))
plot(propagation_grid(i_without_tp)/(2*pi), allfield(:, i_without_tp))
hold on
plot(0:dt:T2mb2, fieldt2mb2)
load HTLSs_harmonic
plot(0:dt:T2mb2, fieldt2mb2)
size(fieldt2mb2)
size(0:dt:T2mb2)
dt
t_ts
dt = 1e-2*2*pi;
save qubits2HOtest_problem fieldw2mb2 fieldt2mb2 u02modes H2modesa H1_2modesa H2_2modesa Hu2modesa H1u2modesa H2u2modesa Hoperations2modesa g2GHz Delta_nuGHz T2mb2
clear all
load qubits2HOtest_problem
Nt_ts = 9;
tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
dt = 1e-2*2*pi;
t_ts = 0.5*(tcheb+1)*dt;
tmidi = round(Nt_ts/2);
test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
% The time-step internal structure including the test-point:
t_ts_tp = [t_ts(1:(Nt_ts - 1)); test_tpoint];
Nt = T2mb2/dt;
dctfactor = T2mb2/(sqrt(Nt*pi));
allt_tp_lasti = Nt*Nt_ts + 1;
allfield = zeros(2, allt_tp_lasti);
for icontrol = 1:2, allfield(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t_ts_tp)/dctfactor; end
propagation_grid = [kron((0:dt:(T2mb2 - dt)), ones(1, Nt_ts)) + kron(ones(1, Nt), [t_ts(1:(Nt_ts - 1)), test_tpoint]), T2mb2];
i_without_tp = 1:allt_tp_lasti;
i_without_tp(Nt_ts:Nt_ts:allt_tp_lasti) = [];
figure
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
for icontrol = 1:2, allfield(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t_ts_tp)/dctfactor; end
propagation_grid = [kron((0:dt:(T2mb2 - dt)), ones(1, Nt_ts)) + kron(ones(1, Nt), [t_ts(1:(Nt_ts - 1)), test_tpoint]), T2mb2];
i_without_tp = 1:allt_tp_lasti;
i_without_tp(Nt_ts:Nt_ts:allt_tp_lasti) = [];
figure
plot(propagation_grid(i_without_tp)/(2*pi), allfield(:, i_without_tp))
hold on
plot(0:dt:T2mb2, fieldt2mb2)
plot((0:dt:T2mb2)/(2*pi), fieldt2mb2)
whos
Hoperations2modesa
[U, mniterc, matvecs, max_errors, history] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 2/3, allfield, [], [-25 35], u02modes, [0, T2mb2], Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
matvecs
Nt
figure
load HTLSs_harmonic
max(max(abs(U-psi2mb2)))
size(psi2mb2)
size(U)
max(abs(U(:, end)-psi2mb2(:,end)))
clear all
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
max_errors
whos
Hu2modesa
nnz(Hu2modesa)
nnz(Hu1modesa)
nnz(H1u1modesa)
nnz(H1u2modesa)
nnz(H2u2modesa)
H1u2modesa
Hoperations2modesa
max(allfield)
max(allfield, 2)
max(allfield, [], 2)
max(abs(allfield), [], 2)
evHu = eig(Hu2modesa)
evH1u = eig(H1u2modesa)
evH2u = eig(H2u2modesa)
30.1+4*1.08
-6.72 - 4*1.08
eig(Hu2modesa - 1.07*H1u2modesa - 1.07*H2u2modesa)
eig(Hu2modesa + 1.07*H1u2modesa + 1.07*H2u2modesa)
[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-10 35], u02modes, [0, T2mb2], Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
max(abs(U(:, end)-U1(:,end)))
max_errors1
max_errors
[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-15 35], u02modes, [0, T2mb2], Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
max_errors1
[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-20 35], u02modes, [0, T2mb2], Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
max_errors1
[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-15 35], u02modes, [0, T2mb2], Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
[Umax, mnitercmax, matvecsmax, max_errorsmax, historymax] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-15 35], u02modes, [0, T2mb2], Nt/2, 9, 9, 1e-6, 10, 16, test_tpoint, false);
max(abs(U(:, end)-Umax(:,end)))
norm(Umax(:, end))
norm(U(:, end))
t01 = 0.5*(tcheb+1)
t01tp = [t01, (t01(tmidi) + t01(tmidi + 1))/2]
for icontrol = 1:2, allfieldmax(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t01tp*T2mb2/140)/(dctfactor*sqrt(2)); end
[Umax, mnitercmax, matvecsmax, max_errorsmax, historymax] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfieldmax, [], [-15 35], u02modes, [0, T2mb2], Nt/2, 9, 9, 1e-6, 10, 16, test_tpoint*2, false);
max(abs(U(:, end)-Umax(:,end)))
max(allfieldmax, [], 2)
max(allfieldmax, [], 2)*2
max(allfield, [], 2)
max(allfieldmax, [], 2)*sqrt(2)
for icontrol = 1:2, allfieldmax(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t01tp*T2mb2/140)/dctfactor; end
[Umax, mnitercmax, matvecsmax, max_errorsmax, historymax] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfieldmax, [], [-15 35], u02modes, [0, T2mb2], Nt/2, 9, 9, 1e-6, 10, 16, test_tpoint*2, false);
max(abs(U(:, end)-Umax(:,end)))
size(allfieldmax)
size(allfield)
Nt
t01tp = [t01(1:(Nt_ts-1)), (t01(tmidi) + t01(tmidi + 1))/2]
t_ts_tp
clear allfieldmax
9*280 + 1
9*140 + 1
allfieldmax = zeros(2, 1261);
for icontrol = 1:2, allfieldmax(icontrol, :) = dctIintgrid1(fieldw2mb2(icontrol, :), T2mb2, t01tp*T2mb2/140)/dctfactor; end
fieldw2mb2
Vt = inv_dctIMtp(fieldw2mb2.', propagation_grid, T2mb2);
max(max(abs(allfield - Vt.')))
max(max(abs(allfield - Vt.'/dctfactor)))
dctfactor
1/dctfactor
figure
plot(propagation_grid(i_without_tp)/(2*pi), Vt(:, i_without_tp).')
plot(propagation_grid(i_without_tp)/(2*pi), Vt(i_without_tp, :).')
figure
plot(propagation_grid(i_without_tp)/(2*pi), Vt(i_without_tp, :).'./allfield)
plot(propagation_grid(i_without_tp)/(2*pi), Vt(i_without_tp, :).'./allfield(:, i_without_tp))
sqrt(N/2)
sqrt(280/2)
sqrt(2/280)
Vt = inv_dctIMtp(fieldw2mb2.', propagation_grid, T2mb2);
Vt = inv_dctIMtp(fieldw2mb2.', propagation_grid, T2mb2, 280);
max(max(abs(allfield - Vt.'/dctfactor)))
tic, for j=1:1e3, Vt = inv_dctIMtp(fieldw2mb2.', propagation_grid, T2mb2, 280); end, toc
tic, for j=1:1e3, allfield_temp = dctIintgrid1(fieldw2mb2(1, :), T2mb2, t_ts_tp)/dctfactor; end, toc
clear Vt allfield_temp
Vw = rand(5, 5);
Vt = rand(5, 5);
Vw = dctIM(Vw);
Vt_back = dctIM(Vw);
Vt - Vt_back
Vw = dctIM(Vt);
Vt_back = dctIM(Vw);
Vt - Vt_back
Vt_back1 = inv_dctIMtp(Vw, 0:4, 4, 4);
Vt_back1-Vt_back
Vt_back1./Vt_back
Vt_back1 = inv_dctIMtp(Vw, 0:4, 4, 4);
Vt_back1-Vt_back
allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:141)*2, T2mb2, 140).'*sqrt(140*pi)/T2mb2;
figure
allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(140*9 + 1))*2, T2mb2, 140).'*sqrt(140*pi)/T2mb2;
140*9+1
2*1260
propagation_grid_max = propagation_grid(1:(140*9 + 1))*2;
propagation_grid_max(end)
propagation_grid_max(end)/(2*pi)
plot(propagation_grid(i_without_tp(1:1261))/(2*pi), allfieldmax(:,i_without_tp(1:1261)))
size(allfieldmax)
plot(propagation_grid(i_without_tp(1:1121))/(2*pi), allfieldmax(:,i_without_tp(1:1121)))
hold on
plot(propagation_grid(i_without_tp)/(2*pi), allfield(:,i_without_tp))
clf
plot(propagation_grid_max(i_without_tp)/(2*pi), allfield(:,i_without_tp))
plot(propagation_grid(i_without_tp)/(2*pi), allfield(:,i_without_tp))
hold on
plot(propagation_grid_max(i_without_tp(1:1121))/(2*pi), allfieldmax(:,i_without_tp(1:1121)))
[Umax, mnitercmax, matvecsmax, max_errorsmax, historymax] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfieldmax, [], [-15 35], u02modes, [0, T2mb2], Nt/2, 9, 9, 1e-6, 10, 16, test_tpoint*2, false);
max(abs(U(:, end)-Umax(:,end)))
max_errorsmax
allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(70*9 + 1))*4, T2mb2, 70).'*sqrt(70*pi)/T2mb2;
figure
save qubits2HOtest_data
factorial(9)^(1/9)
factorial(9)^(1/9)/200
factorial(9)^(1/9)/30
H1u2modesa
H2u2modesa
2.8*2*pi/(factorial(9)^(1/9)/30)
2^9
5^9
3^9
1e5^(1/9)
1e6^(1/9)
2^18
factorial(9)^(1/18)/(40)^9
factorial(9)^(1/18)/(40)^9*1e8^(1/18)
10^8/2^18
10^8/3^18
10^8/2.5^18
factorial(9)^(1/18)/(40)^(1/2)*1e8^(1/18)
factorial(9)^(1/18)/(40)^9*
1e8^(1/18)
2e8^(1/18)
1e8^(1/18)
1e9^(1/18)
factorial(9)/40^9*1e8
(factorial(9)/40^9*1e8)^1/18
(factorial(9)/40^9*1e8)^(1/18)
(factorial(9)/20^9*1e8)^(1/18)
1e-2*2.8*2*pi*35
(factorial(9)/6^9*1e8)^(1/18)
(factorial(9)/6.16^9*1e8)^(1/18)
max_errors
(factorial(9)/6.16^9*2e8)^(1/18)
max_errorsmax
max_errors1
(factorial(9)/6.16^9*5e9)^(1/18)
(factorial(9)/6.16^9/(2e8))^(1/18)
(factorial(9)/6.16^9/(2e-8))^(1/18)
(factorial(9)/6.16^9/(9e-8))^(1/18)
max_errorsmax
allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(140*9 + 1))*2, T2mb2, 140).'*sqrt(140*pi)/T2mb2;
allfieldmax1 = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(70*9 + 1))*4, T2mb2, 70).'*sqrt(70*pi)/T2mb2;
max(abs(U(:, end)-Umax(:,end)))
1 ~= 1:10
[allfieldt, test_tpoint, propagation_grid, propagation_grid_ex_tp] = fieldw2allfieldt(fieldw2mb2, T2mb2, Nt140, 9);
[allfieldt, test_tpoint, propagation_grid, propagation_grid_ex_tp] = fieldw2allfieldt(fieldw2mb2, T2mb2, 140, 9);
max(abs(allfieldt-allfield))
size(allfieldt)
max(abs(allfieldt-allfieldmax))
max(abs(allfieldt-allfieldmax), [], 2)
test_tpoint
test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
test_tpoint*2
propagation_grid = [kron((0:dt:(T2mb2 - dt)), ones(1, Nt_ts)) + kron(ones(1, Nt), [t_ts(1:(Nt_ts - 1)), test_tpoint]), T2mb2];
[allfieldtemp, test_tpointtemp, propagation_grid_temp, propagation_grid_ex_tp_temp] = fieldw2allfieldt(fieldw2mb2, T2mb2, 140, 9);
propagation_grid_temp(1:30)./propagation_grid(1:30)
propagation_grid_ex_tp_temp
propagation_grid_ex_tp_temp(1:30)
figure
plot(propagation_grid_max_ex_tp/(2*pi), allfieldmax(:,i_without_tp(1:1121)))
plot(propagation_grid_ex_tp_temp/(2*pi), allfieldmax(:,i_without_tp(1:1121)))
[allfieldtemp, test_tpointtemp, propagation_grid_temp, propagation_grid_ex_tp_temp, i_ex_tp_temp] = fieldw2allfieldt(fieldw2mb2, T2mb2, 140, 9);
plot(propagation_grid_ex_tp_temp/(2*pi), allfieldmax(:,i_ex_tp_temp))
clear allfieldtemp test_tpointtemp propagation_grid_temp propagation_grid_ex_tp_temp i_ex_tp_temp
[allfieldt70, test_tpoint70, propagation_grid70, propagation_grid70, i_ex_tp70] = fieldw2allfieldt(fieldw2mb2, T2mb2, 70, 9);
max(abs(allfieldt70-allfieldmax1), [], 2)
clear allfieldmax1
[U70, mniterc70, matvecs70, max_errors70, history70] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfieldt70, [], [-15 35], u02modes, [0, T2mb2], 70, 9, 9, 1e-6, 10, 16, test_tpoint70, false);
vecnorm(U70)
max_errors70
(factorial(9)/6.16^9/(3e-6))^(1/18)
(factorial(9)/6.16^9/(9e-8))^(1/18)
4e-2*2.8*2*pi*35
2.8*2*pi/70
2.8*2*pi/70*35
T2mb2/70*35
2.8*2*pi/70
ans/(2*pi)
4e-2*2.8*2*pi
2.8*2*pi
ans/100
(factorial(9)/8.7965^9/(3e-6))^(1/18)
7/5
[allfield50, test_tpoint70, propagation_grid70, propagation_grid70, i_ex_tp70] = fieldw2allfieldt(fieldw2mb2, T2mb2, 50, 9);
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, [0, T2mb2], 70, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
[allfieldt70, test_tpoint70, propagation_grid70, propagation_grid70, i_ex_tp70] = fieldw2allfieldt(fieldw2mb2, T2mb2, 70, 9);
[allfield50, test_tpoint50, propagation_grid50, propagation_grid50, i_ex_tp50] = fieldw2allfieldt(fieldw2mb2, T2mb2, 50, 9);
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, [0, T2mb2], 70, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
size(allfield50)
test_tpoint50
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, [0, T2mb2], 70, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, [0, T2mb2], 50, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
max_errors
max_errors50
(factorial(9)/4.96^9/(4.96e-6))^(1/18
2.8*2*pi/50*35
(factorial(9)/12.3^9/(4.96e-6))^(1/18
(factorial(9)/12.3^9/(4.96e-6))^(1/18)
(factorial(9)/12.3^9/(4.96e-5))^(1/18)
[allfield25, test_tpoint25, propagation_grid25, propagation_grid25, i_ex_tp25] = fieldw2allfieldt(fieldw2mb2, T2mb2, 25, 9);
[U25, mniterc25, matvecs25, max_errors25, history25] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield25, [], [-15 35], u02modes, [0, T2mb2], 25, 9, 9, 1e-6, 10, 16, test_tpoint25, false);
vecnorm(U250)
vecnorm(U25)
[U25, mniterc25, matvecs25, max_errors25, history25] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield25, [], [-15 35], u02modes, propagation_grid25, 25, 9, 9, 1e-6, 10, 16, test_tpoint25, false);
[P, D] = eig(Hu2modesa);
eig(Hu2modesa);
eig(Hu2modesa)
[allfield40, test_tpoint40, propagation_grid40, propagation_grid40, i_ex_tp40] = fieldw2allfieldt(fieldw2mb2, T2mb2, 40, 9);
[U40, mniterc40, matvecs40, max_errors40, history40] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield40, [], [-15 35], u02modes, propagation_grid40, 40, 9, 9, 1e-6, 10, 16, test_tpoint40, false);
vecnorm(U40)
[allfield45, test_tpoint45, propagation_grid45, propagation_grid45, i_ex_tp45] = fieldw2allfieldt(fieldw2mb2, T2mb2, 45, 9);
[U45, mniterc45, matvecs45, max_errors45, history45] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield45, [], [-15 35], u02modes, propagation_grid45, 45, 9, 9, 1e-6, 10, 16, test_tpoint45, false);
vecnorm(U45)
max_errors45
2.8*2*pi/45*35
(factorial(9)/13.68^9/(1.176e-4))^(1/18)
save qubits2HOtest_data
whos
%-- 21/07/2022 23:28 --%
load qubits2HOtest_problem
load qubits2HOtest_data
(factorial(9)/10^9/(1e-7))^(1/18)
mniter70
mniterc70
mniterc50
mnitercmax
Nt
mniterc
H1u2modesa
H2u2modesa
2.8*2*pi/50
(2.8*2*pi/50)^2
(2.8*2*pi/50)^2*2
%-- 24/07/2022 10:00 --%
load qubits2HOtest_data
doc sparse
[P, D] = eig(full(Hu2modesa));
E = diag(D)
U50E = P\U50;
figure
plot(0:T2mb2/50:T2mb2, U50E)
size(U50E)
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
size(U50)
size(propgation_grid50)
size(propagation_grid50)
[allfield50, test_tpoint50, propagation_grid50, propagation_grid50ex_tp, i_ex_tp50] = fieldw2allfieldt(fieldw2mb2, T2mb2, 50, 9);
[allfield45, test_tpoint45, propagation_grid45, propagation_grid45ex_tp, i_ex_tp45] = fieldw2allfieldt(fieldw2mb2, T2mb2, 45, 9);
[allfield40, test_tpoint40, propagation_grid40, propagation_grid40ex_tp, i_ex_tp40] = fieldw2allfieldt(fieldw2mb2, T2mb2, 40, 9);
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
mniterc
mniterc50
size(propagation_grid50)
size(U50)
U50E = P\U50;
plot(0:T2mb2/50:T2mb2, U50E)
[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-15 35], u02modes, 0:T2mb2/280:T2mb2, Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
Nt
U1E = P\U1;
figure
plot(0:T2mb2/280:T2mb2, U1E)
norm(U1(:,end)-U50(:,end))
[U45, mniterc45, matvecs45, max_errors45, history45] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield45, [], [-15 35], u02modes, propagation_grid45(1:9:end), 45, 9, 9, 1e-6, 10, 16, test_tpoint45, false);
norm(U1(:,end)-U45(:,end))
U45E = P\U45;
figure
plot(0:T2mb2/45:T2mb2, U45E)
[U40, mniterc40, matvecs40, max_errors40, history40] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield40, [], [-15 35], u02modes, propagation_grid40(1:9:end), 40, 9, 9, 1e-6, 10, 16, test_tpoint40, false);
U40E = P\U40;
figure
plot(0:T2mb2/40:T2mb2, U40E.conj(U40E))
plot(0:T2mb2/40:T2mb2, U40E.*conj(U40E))
figure
plot(0:T2mb2/280:T2mb2, log10(U1E.*conj(U1E)))
figure
plot(0:T2mb2/45:T2mb2, log10(U45E.*conj(U45E)))
figure
plot(0:T2mb2/50:T2mb2, log10(U50E.*conj(U50E)))
figure
plot(0:T2mb2/40:T2mb2, log10(U40E.*conj(U40E)))
[E, U40E(:, end).*conj(U40E(:, end))]
U40(:, end).*conj(U40(:, end))
fieldw2mb2
dw
dw = T2mb2/(2*pi*2.8)
dw = pi/(2*pi*2.8)
wgrid = 0:dw:215*dw;
allfield40temp = inv_dctIMtp_wgrid(fieldw2mb2, propagation_grid40, wgrid, 280);
allfield40temp = inv_dctIMtp_wgrid(fieldw2mb2.', propagation_grid40, wgrid, 280);
allfield40temp = inv_dctIMtp_wgrid(fieldw2mb2(:, 1:216).', propagation_grid40, wgrid, 280);
max(abs(allfield40 - allfield40temp))
size(allfield40temp)
size(allfield40)
allfield40temp = inv_dctIMtp_wgrid(fieldw2mb2(:, 1:216).', propagation_grid40, wgrid, 280).';
max(abs(allfield40 - allfield40temp))
max(abs(allfield40 - allfield40temp), [], 2)
allfield40temp = inv_dctIMtp_wgrid(fieldw2mb2(:, 1:216).', propagation_grid40, wgrid, 280).'*sqrt;
allfield40temp = inv_dctIMtp_wgrid(fieldw2mb2(:, 1:216).', propagation_grid40, wgrid, 280).'*sqrt(280*pi)/(2*pi*2.8);
max(abs(allfield40 - allfield40temp), [], 2)
Gop2mb2 = @(u, t, v) -1i*(Hu2modesa*v - inv_dctIMtp_wgrid(fieldw2mb2(1, 1:216).', t, wgrid, 280).'/dctfactor*(H1u2modesa*v) - inv_dctIMtp_wgrid(fieldw2mb2(2, 1:216).', t, wgrid, 280).'/dctfactor*(H2u2modesa*v))
dctfactor
T2mb2/(sqrt(pi*280))
Gdiff_op2mb2 = @(u1, t1, u2, t2) 1i*(((inv_dctIMtp_wgrid(fieldw2mb2(1,1:216).',t1,wgrid,280).' - inv_dctIMtp_wgrid(fieldw2mb2(1,1:216).',t2,wgrid,280).')/dctfactor).*(H1u2modesa*u1) + ((inv_dctIMtp_wgrid(fieldw2mb2(2,1:216).',t1,wgrid,280).' - inv_dctIMtp_wgrid(fieldw2mb2(2,1:216).',t2,wgrid,280).')/dctfactor).*(H2u2modesa*u1))
[U50temp, mniterc50temp, matvecs50temp, max_errors50temp, history50temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
[U50temp, mniterc50temp, matvecs50temp, max_errors50temp, history50temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, false);
Gop2mb2 = @(u, t, v) -1i*(Hu2modesa*v - [H1u2modesa*v, H2u2modesa*v]*(inv_dctIMtp_wgrid(fieldw2mb2(:, 1:216).', t, wgrid, 280).'/dctfactor))
fieldw2mb2T = fieldw2mb2(:,1:216).';
[U50temp, mniterc50temp, matvecs50temp, max_errors50temp, history50temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, false);
mniterc50temp
Gdiff_op2mb2 = @(u1, t1, u2, t2) 1i*(((inv_dctIMtp_wgrid(fieldw2mb2T(:, 1),t1,wgrid,280).' - inv_dctIMtp_wgrid(fieldw2mb2T(:, 1).',t2,wgrid,280))/dctfactor).*(H1u2modesa*u1) + ((inv_dctIMtp_wgrid(fieldw2mb2T(:, 2),t1,wgrid,280).' - inv_dctIMtp_wgrid(fieldw2mb2T(:, 2),t2,wgrid,280))/dctfactor).*(H2u2modesa*u1))
Gop2mb2 = @(u, t, v) -1i*(Hu2modesa*v - [H1u2modesa*v, H2u2modesa*v]*(inv_dctIMtp_wgrid(fieldw2mb2T, t, wgrid, 280).'/dctfactor))
[U50temp, mniterc50temp, matvecs50temp, max_errors50temp, history50temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, false);
Gdiff_op2mb2 = @(u1, t1, u2, t2) 1i*(((inv_dctIMtp_wgrid(fieldw2mb2T(:, 1),t1,wgrid,280).' - inv_dctIMtp_wgrid(fieldw2mb2T(:, 1),t2,wgrid,280))/dctfactor).*(H1u2modesa*u1) + ((inv_dctIMtp_wgrid(fieldw2mb2T(:, 2),t1,wgrid,280).' - inv_dctIMtp_wgrid(fieldw2mb2T(:, 2),t2,wgrid,280))/dctfactor).*(H2u2modesa*u1))
[U50temp, mniterc50temp, matvecs50temp, max_errors50temp, history50temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, false);
[U1temp, mniterc1temp, matvecs1temp, max_errors1temp, history1temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid(1:9:end), 280, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U1temp-U1)))
tic,[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-15 35], u02modes, 0:T2mb2/280:T2mb2, Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);toc
tic,[U1temp, mniterc1temp, matvecs1temp, max_errors1temp, history1temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid(1:9:end), 280, 9, 9, 1e-6, 10, 16, false);toc
max_errors1temp
max_errors1
[U70temp, mniterc70temp, matvecs70temp, max_errors70temp, history70temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid70(1:9:end), 70, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U70temp-U70)))
(max(abs(U70temp(:, end)-U70(:,end))))
allfield40temp1 = inv_dctIMtp_wgrid(fieldw2mb2T, propagation_grid40, wgrid, 280).'/dctfactor;
max(abs(allfield40temp1 - allfield40temp), [], 2)
get_fieldt = @(t) inv_dctIMtp_wgrid(fieldw2mb2T, t, wgrid, 280).'/dctfactor;
get_fieldt1 = @(t) inv_dctIMtp_wgrid(fieldw2mb2T(:,1), t, wgrid, 280).'/dctfactor;
get_fieldt2 = @(t) inv_dctIMtp_wgrid(fieldw2mb2T(:,2), t, wgrid, 280).'/dctfactor;
tic, allfieldmax = inv_dctIMtp(fieldw2mb2.', propagation_grid(1:(140*9 + 1))*2, T2mb2, 140).'*sqrt(140*pi)/T2mb2;toc
tic, allfieldtemp = inv_dctIMtp_wgrid(fieldw2mb2T, propagation_grid, T2mb2, 280).'/dctfactor;toc
tic, allfieldtemp = inv_dctIMtp_wgrid(fieldw2mb2T, propagation_grid, wgrid, 280).'/dctfactor;toc
max(abs(get_fieldt(propagation_grid40) - allfield40temp), [], 2)
Gop2mb2 = @(u, t, v) -1i*(Hu2modesa*v - [H1u2modesa*v, H2u2modesa*v]*get_fieldt(t))
Gdiff_op2mb2 = @(u1, t1, u2, t2) 1i*((get_fieldt1(t1) - get_fieldt1(t2)).*(H1u2modesa*u1) + (get_fieldt2(t1) - get_fieldt2(t2)).*(H2u2modesa*u1))
[U1temp, mniterc1temp, matvecs1temp, max_errors1temp, history1temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-15 35], u02modes, propagation_grid(1:9:end), 280, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U1temp-U1)))
[U1temp, mniterc1temp, matvecs1temp, max_errors1temp, history1temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid(1:9:end), 280, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U1temp-U1)))
[U1temp, mniterc1temp, matvecs1temp, max_errors1temp, history1temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid(1:9:end), 280, 9, 9, 1e-6, 10, 16, false);
tic, [U1temp, mniterc1temp, matvecs1temp, max_errors1temp, history1temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid(1:9:end), 280, 9, 9, 1e-6, 10, 16, false);toc
[U50temp, mniterc50temp, matvecs50temp, max_errors50temp, history50temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U50temp-U50)))
[U70temp, mniterc70temp, matvecs70temp, max_errors70temp, history70temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid70(1:9:end), 70, 9, 9, 1e-6, 10, 16, false);
max_errors70
max_errors70temp
max(max(abs(U70temp-U70)))
max(max(abs(U70temp(:,end)-U70(:,end))))
max(max(abs(U70temp(:,end)-U1(:,end))))
max(max(abs(U70(:,end)-U1(:,end))))
max(max(abs(U50(:,end)-U1(:,end))))
max(max(abs(U50temp(:,end)-U1(:,end))))
[U70temp, mniterc70temp, matvecs70temp, max_errors70temp, history70temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid70(1:9:end), 70, 9, 9, 1e-6, 10, 16, false);
[allfield70, test_tpoint70, propagation_grid70, propagation_grid70ex_tp, i_ex_tp70] = fieldw2allfieldt(fieldw2mb2, T2mb2, 70, 9);
[U70temp, mniterc70temp, matvecs70temp, max_errors70temp, history70temp] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, propagation_grid70(1:9:end), 70, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U70temp(:,end)-U1(:,end))))
max(max(abs(U70temp(:,end)-U70(:,end))))
[allfield50, test_tpoint50, propagation_grid50, propagation_grid50ex_tp, i_ex_tp50] = fieldw2allfieldt(fieldw2mb2, T2mb2, 50, 9);
[U50, mniterc50, matvecs50, max_errors50, history50] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, 1e-6, 10, 16, test_tpoint50, false);
max(max(abs(U50temp(:,end)-U50(:,end))))
max(history50.niter - history50temp.niter)
max(history70.niter - history70temp.niter)
max_errors70temp
max_errors70
[U43, mniterc43, matvecs43, max_errors43, history43] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/43:T2mb2, 43, 9, 9, 1e-6, 10, 16, false);
max(max(abs(U43(:,end)-U1(:,end))))
U43E = P\U43;
figure
plot(0:T2mb2/43:T2mb2, log10(U43E.*conj(U43E)))
[U42, mniterc42, matvecs42, max_errors42, history42] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/42:T2mb2, 42, 9, 9, 1e-6, 10, 16, false);
figure
plot(0:T2mb2/42:T2mb2, log10(U42E.*conj(U42E)))
U42E = P\U42;
plot(0:T2mb2/42:T2mb2, log10(U42E.*conj(U42E)))
max_errors43
(35*1.6957e-4)
(factorial(9)*43^9/(T2mb2*35)^9/(1.176e-4))^(1/18)
(factorial(9)*43^9/(T2mb2*35)^9/(1.7e-4))^(1/18)
43/ans
max_errors1
(factorial(9)*(280/(T2mb2*35))^9/(9.21e-8))^(1/18)
(factorial(9)*(43/(T2mb2*35))^9/(1.7e-4))^(1/18)
(factorial(9)*(280/(T2mb2*35))^9/(9.21e-8))^(1/18)
280/ans
max_errorsmax
(factorial(9)*(140/(T2mb2*35))^9/(6.326e-9))^(1/18)
140/ans
[U1, mniterc1, matvecs1, max_errors1, history1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield, [], [-15 35], u02modes, 0:T2mb2/280:T2mb2, Nt, 9, 9, 1e-6, 10, 16, test_tpoint, false);
max_errors1
whos
clear U1temp mniterc1temp matvecs1temp max_errors1temp history1temp
clear U50temp mniterc50temp matvecs50temp max_errors50temp history50temp
clear U70temp mniterc70temp matvecs70temp max_errors70temp history70temp
whos
save qubits2HOtest_data
[U50i, mniterc50i, matvecs50i, max_errors50i, history50i] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield50, [], [-15 35], u02modes, propagation_grid50(1:9:end), 50, 9, 9, eps, 1, 16, [], false);
sqnorm(U50i(:,end))
[U70i, mniterc70i, matvecs70i, max_errors70i, history70i] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfieldt70, [], [-15 35], u02modes, [0, T2mb2], 70, 9, 9, eps, 1, 16, test_tpoint70, false);
mniterc50i
mniterc70i
sqnorm(U70i(:,end))
[U80i, mniterc80i, matvecs80i, max_errors80i, history80i] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/80:T2mb2, 80, 9, 9, eps, 1, 16, false);
mniterc80i
sqnorm(U80i(:,end))
[U75i, mniterc75i, matvecs75i, max_errors75i, history75i] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/80:T2mb2, 75, 9, 9, eps, 1, 16, false);
sqnorm(U75i(:,end))
U80iE = P\U80i;
figure
plot(0:T2mb2/80:T2mb2, log10(U80E.*conj(U80E)))
plot(0:T2mb2/80:T2mb2, log10(U80iE.*conj(U80iE)))
U75iE = P\U75i;
plot(0:T2mb2/75:T2mb2, log10(U75iE.*conj(U75iE)))
size(U75i)
[U75i, mniterc75i, matvecs75i, max_errors75i, history75i] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/75:T2mb2, 75, 9, 9, eps, 1, 16, false);
plot(0:T2mb2/75:T2mb2, log10(U75iE.*conj(U75iE)))
U75iE = P\U75i;
plot(0:T2mb2/75:T2mb2, log10(U75iE.*conj(U75iE)))
max_errors
plot((0:T2mb2/40:T2mb2)/(2*pi), log10(U40E.*conj(U40E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\left|\left<\phi_n|\psi(T)\right>\right|$', 'interpreter', 'latex')
ylabel('$\left|\left<\varphi_n|\psi(T)\right>\right|^2$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2$\right)', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
title('\Delta t = T/40', 'interpreter', 'latex')
title('$\Delta t = T/40$', 'interpreter', 'latex')
plot((0:T2mb2/280:T2mb2)/(2*pi), log10(U1E.*conj(U1E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
title('$\Delta t = T/280$', 'interpreter', 'latex')
figure
plot((0:T2mb2/43:T2mb2)/(2*pi), log10(U43E.*conj(U43E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
title('$\Delta t = T/43$', 'interpreter', 'latex')
plot((0:T2mb2/42:T2mb2)/(2*pi), log10(U42E.*conj(U42E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
title('$\Delta t = T/42$', 'interpreter', 'latex')
save qubits2HOtest_data
[allfield500, test_tpoint500, propagation_grid500, propagation_grid500ex_tp, i_ex_tp500] = fieldw2allfieldt(fieldw2mb2, T2mb2, 500, 9);
[U500, mniterc500, matvecs500, max_errors500, history500] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield500, [], [-15 35], u02modes, propagation_grid500(1:9:end), 500, 9, 9, eps, 10, 16, [], false);
max_errors500
[U500_1, mniterc500_1, matvecs500_1, max_errors500_1, history500_1] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield500, [], [-15 35], u02modes, propagation_grid500(1:9:end), 500, 9, 13, eps, 10, 16, [], false);
max_errors500_1
max(abs(U500(:,end)-U500_1(:,end)))
clear history500 history500_1
[allfield400, test_tpoint400, propagation_grid400, propagation_grid400ex_tp, i_ex_tp400] = fieldw2allfieldt(fieldw2mb2, T2mb2, 400, 9);
[U400, mniterc400, matvecs400, max_errors400, history400] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield400, [], [-15 35], u02modes, propagation_grid400(1:9:end), 400, 9, 13, eps, 10, 16, [], false);
max(abs(U400(:,end)-U500_1(:,end)))
max(abs(U1(:,end)-U500_1(:,end)))
Uex = U500_1(:,end);
clear history400
save Uex_2q2rtest Uex fieldw2mb2 dctfactor T2mb2
log10(400)
log10(80)
log10(500)
T2mb2/sqrt(280*pi)
save Uex_2q2rtest Uex fieldw2mb2T
save Uex_2q2rtest Uex fieldw2mb2
save Uex_2q2rtest Uex
[allNt9, allmv9, aller9, max_ers9] = error2q2r(T2mb2, 9, 9, 80, 8);
[allNt9, allmv9, aller9, max_ers9] = error2q2r(T2mb2, 9, 9, 80, 14);
polyfit(log10(allmv9), log10(aller), 1)
polyfit(log10(allmv9), log10(aller9), 1)
polyfit(log10(allmv9(1:4)), log10(aller9(1:4)), 1)
polyfit(log10(allmv9(5:11)), log10(aller9(5:11)), 1)
log10(500)
allNt9
[U80i7, mniterc80i7, matvecs80i7, max_errors80i7, history80i7] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/80:T2mb2, 80, 7, 7, eps, 1, 16, false);
sqnorm(U80i7(:,end))
sqnorm(U80i(:,end))
max_errors80i
max_errors80i7
(factorial(7)*(80/(T2mb2*35))^7/(1.31e-4))^(1/14)
80/ans
[U64i7, mniterc64i7, matvecs64i7, max_errors64i7, history64i7] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/64:T2mb2, 64, 7, 7, eps, 1, 16, false);
sqnorm(U64i7(:,end))
[U60i7, mniterc60i7, matvecs60i7, max_errors60i7, history60i7] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/60:T2mb2, 60, 7, 7, eps, 1, 16, false);
sqnorm(U60i7(:,end))
[allNt9, allmv9, aller9, max_ers9] = error2q2r(T2mb2, 9, 9, 80, 14, 1, True);
[allNt9, allmv9, aller9, max_ers9] = error2q2r(T2mb2, 9, 9, 80, 14, 1, true);
figure
[allNt9, allmv9, aller9, max_ers9] = error2q2r(T2mb2, 9, 9, 80, 14);
[allNt9ar, allmv9ar, aller9ar, max_ers9ar] = error2q2r(T2mb2, 9, 9, 80, 14, 1, true);
plot(log10(allmv9), log10(aller9), '-o')
hold on
plot(log10(allmv9ar), log10(aller9ar), '-o')
[allNt7, allmv7, aller7, max_ers7] = error2q2r(T2mb2, 7, 7, 64, 14);
[allNt7, allmv7, aller7, max_ers7] = error2q2r(T2mb2, 7, 7, 64, 16);
[allNt7, allmv7, aller7, max_ers7] = error2q2r(T2mb2, 7, 7, 64, 17);
polyfit(log10(allmv9(1:4)), log10(aller9(1:4)), 1)
polyfit(log10(allmv7(4:11)), log10(aller7(4:11)), 1)
polyfit(log10(allmv7(1:4)), log10(aller7(1:4)), 1)
max_errors64i7
figure
plot(log10(allmv9), log10(aller9), '-o')
hold on
plot(log10(allmv7), log10(aller7), '-o')
[U64i97, mniterc64i97, matvecs64i97, max_errors64i97, history64i97] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/64:T2mb2, 64, 9, 7, eps, 1, 16, false);
sqnorm(U64i97(:,end))
[U75i97, mniterc75i97, matvecs75i97, max_errors75i97, history75i97] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/75:T2mb2, 75, 9, 7, eps, 1, 16, false);
sqnorm(U75i97(:,end))
[allNt97, allmv97, aller97, max_ers97] = error2q2r(T2mb2, 9, 7, 80, 17);
plot(log10(allmv97), log10(aller97), '-o')
[allNt95, allmv95, aller95, max_ers95] = error2q2r(T2mb2, 9, 5, 80, 17);
plot(log10(allmv95), log10(aller95), '-o')
[allNt93, allmv93, aller93, max_ers93] = error2q2r(T2mb2, 9, 3, 80, 17);
plot(log10(allmv93(2:end)), log10(aller93(2:end)), '-o')
[allNt75, allmv75, aller75, max_ers75] = error2q2r(T2mb2, 7, 5, 64, 17);
plot(log10(allmv75(3:end)), log10(aller75(3:end)), '-o')
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/500);
sqnorm(URK(:,end))
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/100);
sqnorm(URK(:,end))
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/300);
sqnorm(URK(:,end))
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/200);
sqnorm(URK(:,end))
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/150);
sqnorm(URK(:,end))
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/180);
sqnorm(URK(:,end))
URK = RK4uf(@(t, u) Gop2mb2(u, t, u), [0 T2mb2], u02modes, T2mb2/200);
URK = RK4uf(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/200);
sqnorm(URK(:,end))
dw
[allNtRK, allmvRK, allerRK] =  error2q2rRK(T2mb2, 200, 20);
allmvRK
hold on
plot(log10(allmvRK), log10(allerRK), '-o')
figure
plot(log10(allmvRK), log10(allerRK), '-o')
hold on
plot(log10(allmv75(3:end)), log10(aller75(3:end)), '-o')
[U64i55, mniterc64i55, matvecs64i5, max_errors64i5, history64i5] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/64:T2mb2, 64, 5, 5, eps, 1, 16, false);
sqnorm(U64i55(:,end))
(factorial(5)*(64/(T2mb2*35))^5/(3.7e-2))^(1/10)
64/ans
[U90i55, mniterc90i55, matvecs90i5, max_errors90i5, history90i5] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/90:T2mb2, 90, 5, 5, eps, 1, 16, false);
sqnorm(U90i55(:,end))
(factorial(5)*(64/(T2mb2*35))^5/(6.15e-3))^(1/10)
(factorial(5)*(90/(T2mb2*35))^5/(6.15e-3))^(1/10)
[U95i55, mniterc95i55, matvecs95i5, max_errors95i5, history95i5] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/95:T2mb2, 95, 5, 5, eps, 1, 16, false);
sqnorm(U95i55(:,end))
save qubits2HOtest_data
[U90_55, mniterc90_55, matvecs90_5, max_errors90_5, history90_5] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/90:T2mb2, 90, 5, 5, 1e-6, 10, 16, false);
sqnorm(U90_55(:,end))
clear U90_55 mniterc90_55 matvecs90_5 max_errors90_5 history90_5
[U100i5, mniterc100i5, matvecs100i5, max_errors100i5, history100i5] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/100:T2mb2, 100, 5, 5, eps, 1, 16, false);
sqnorm(U100i5(:,end))
[allNtRK, allmvRK, allerRK] =  error2q2rRK(T2mb2, 200, 40);
figure
plot(log10(allmvRK(1:28)), log10(allerRK(1:28)), '-o')
hold on
plot(log10(allmv97), log10(aller97), '-o')
plot(log10(allmv97(1:13)), log10(aller97(1:13)), '-o')
plot(log10(allmv7), log10(aller7), '-o')
plot(log10(allmv95), log10(aller95), '-o')
plot(log10(allmv93), log10(aller93), '-o')
plot(log10(allmv95), log10(aller95), '-o')
plot(log10(allmv95(1:13)), log10(aller95(1:13)), '-o')
plot(log10(allmv93(2:end)), log10(aller93(2:end)), '-o')
figure
plot(log10(allmv97(1:13)), log10(aller97(1:13)), '-o')
hold on
plot(log10(allmv93(2:end)), log10(aller93(2:end)), '-o')
polyfit(log10(allmv93(1:11)), log10(aller93(1:11)), 1)
polyfit(log10(allmv93(2:12)), log10(aller93(2:12)), 1)
[allNt94, allmv94, aller94, max_ers94] = error2q2r(T2mb2, 9, 4, 80, 17);
polyfit(log10(allmv94(2:11)), log10(aller94(2:11)), 1)
hold on
plot(log10(allmv95(1:13)), log10(aller95(1:13)), '-o')
plot(log10(allmv93(2:end)), log10(aller93(2:end)), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
polyfit(log10(allmvRK(5:27)), log10(aller94(5:27)), 1)
polyfit(log10(allmvRK(5:27)), log10(allerRK(5:27)), 1)
save qubits2HOtest_data
[allNt95, allmv95, aller95, max_ers95] = error2q2r(T2mb2, 9, 5, 80, 17, 1, true);
hold on
plot(log10(allmv95(1:13)), log10(aller95(1:13)), '-o')
[allNt95, allmv95, aller95, max_ers95] = error2q2r(T2mb2, 9, 5, 80, 17);
[allNt95ar, allmv95ar, aller95ar, max_ers95ar] = error2q2r(T2mb2, 9, 5, 80, 17, 1, true);
hold on
plot(log10(allmv95(1:13)), log10(aller95(1:13)), '-o')
[allNt93ar, allmv93ar, aller93ar, max_ers93ar] = error2q2r(T2mb2, 9, 3, 80, 17, 1, true);
hold on
plot(log10(allmv93(2:end)), log10(aller93(2:end)), '-o')
hold on
plot(log10(allmv93ar(1:13)), log10(aller93ar(1:13)), '-o')
10^0.5
10^0.7
10^0.8
10^0.9
HopPWC = @(v, fields) Hu2modesa*v - [H1u2modesa*v, H2u2modesa*v]*fields;
tgridPWC200 = (T2mb2/200/2:T2mb2/200:(T2mb2 - T2mb2/200/2)).';
fieldPWC200 = inv_dctIMtp_wgrid(fieldw2mb2T, tgridPWC200, wgrid, 280).'/dctfactor;
T2mb2/200*50
UPWC200 = SchrPWCcheb(Hop, u02modes, fieldPWC200, [-15 35], T2mb2, 200, 9);
UPWC200 = SchrPWCcheb(HopPWC, u02modes, fieldPWC200, [-15 35], T2mb2, 200, 9);
sqnorm(UPWC200(:,end))
UPWC200_2 = SchrPWCcheb(HopPWC, u02modes, fieldPWC200, [-15 35], T2mb2, 200, 4);
sqnorm(UPWC200_2(:,end))
UPWC200_2 = SchrPWCcheb(HopPWC, u02modes, fieldPWC200, [-15 35], T2mb2, 200, 7);
sqnorm(UPWC200_2(:,end))
UPWC200_2 = SchrPWCcheb(HopPWC, u02modes, fieldPWC200, [-15 35], T2mb2, 200, 11);
sqnorm(UPWC200_2(:,end))
UPWC200_2 = SchrPWCcheb(HopPWC, u02modes, fieldPWC200, [-15 35], T2mb2, 200, 13);
sqnorm(UPWC200_2(:,end))
4-sqnorm(UPWC200_2(:,end))
norm(Uex(:,end)-UPWC200_2(:,end))
norm(Uex(:,end)-UPWC200(:,end))
tgridPWC100 = (T2mb2/100/2:T2mb2/100:(T2mb2 - T2mb2/100/2)).';
UPWC100 = SchrPWCcheb(HopPWC, u02modes, fieldPWC100, [-15 35], T2mb2, 100, 13);
fieldPWC100 = inv_dctIMtp_wgrid(fieldw2mb2T, tgridPWC100, wgrid, 280).'/dctfactor;
UPWC100 = SchrPWCcheb(HopPWC, u02modes, fieldPWC100, [-15 35], T2mb2, 100, 13);
sqnorm(UPWC100(:,end))
UPWC100 = SchrPWCcheb(HopPWC, u02modes, fieldPWC100, [-15 35], T2mb2, 100, 17);
sqnorm(UPWC100(:,end))
4-sqnorm(UPWC100(:,end))
norm(Uex(:,end)-UPWC100(:,end))
aller7
tgridPWC50 = (T2mb2/50/2:T2mb2/50:(T2mb2 - T2mb2/50/2)).';
fieldPWC50 = inv_dctIMtp_wgrid(fieldw2mb2T, tgridPWC50, wgrid, 280).'/dctfactor;
size(UPWC100)
UPWC50 = SchrPWCcheb(HopPWC, u02modes, fieldPWC50, [-15 35], T2mb2, 100, 20);
UPWC50 = SchrPWCcheb(HopPWC, u02modes, fieldPWC50, [-15 35], T2mb2, 50, 20);
sqnorm(UPWC50(:,end))
UPWC50 = SchrPWCcheb(HopPWC, u02modes, fieldPWC50, [-15 35], T2mb2, 50, 25);
sqnorm(UPWC50(:,end))
4-sqnorm(UPWC50(:,end))
norm(Uex(:,end)-UPWC50(:,end))
UPWC50 = SchrPWCcheb(HopPWC, u02modes, fieldPWC50, [-15 35], T2mb2, 50, 15);
norm(Uex(:,end)-UPWC50(:,end))
UPWC50 = SchrPWCcheb(HopPWC, u02modes, fieldPWC50, [-15 35], T2mb2, 50, 20);
norm(Uex(:,end)-UPWC50(:,end))
T2mb2/50*50/2
[allNtPWC20, allmvPWC20, allerPWC20] =  error2q2rPWC(T2mb2, 20, 50, 30);
polyfit(log10(allmvPWC20), log10(allerPWC20), 1)
[allNtPWC13, allmvPWC13, allerPWC13] =  error2q2rPWC(T2mb2, 13, 100, 30);
polyfit(log10(allmvPWC13), log10(allerPWC13), 1)
norm(Uex(:,end)-UPWC200(:,end))
norm(Uex(:,end)-UPWC200_2(:,end))
[allNtPWC10, allmvPWC10, allerPWC10] =  error2q2rPWC(T2mb2, 13, 200, 30);
polyfit(log10(allmvPWC10), log10(allerPWC10), 1)
[allNtPWC10, allmvPWC10, allerPWC10] =  error2q2rPWC(T2mb2, 10, 200, 30);
polyfit(log10(allmvPWC10), log10(allerPWC10), 1)
figure
plot(log10(allmvRK(1:28)), log10(allerRK(1:28)), '-o')
hold on
plot(log10(allmvPWC10), log10(allerPWC10), '-o')
plot(log10(allmvPWC13), log10(allerPWC13), '-o')
plot(log10(allmvPWC20), log10(allerPWC20), '-o')
plot(log10(allmv7), log10(aller7), '-o')
[allNtPWC20Hb, allmvPWC20Hb, allerPWC20Hb] =  error2q2rPWC_Hb(T2mb2, 20, 50, 30);
hold on
plot(log10(allmvPWC20), log10(allerPWC20), '-o')
polyfit(log10(allmvPWC20Hb), log10(allerPWC20Hb), 1)
save qubits2HOtest_data
xlabel('log(matvecs)')
ylabel('log(error)')
0.2/g2GHz
g2GHz
0.25/0.03
sum([41, 0, 0, 216, 27, 272, 27, 216, 41])
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/200);
ti
fderiv(t, U(:, ti))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/200);
sqnorm(URK7(:,end))
sqnorm(URK7(:,2))
sqnorm(URK7(:,3))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/500);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/1e3);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/2e3);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/5e3);
sqnorm(URK7(:,end))
norm(Uex(:,end)-URK7(:,end))
70/5
255/5
265/5
379/7
849/3
833/7
156/6
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/200);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/500);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/5000);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/5000);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/5000);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/200);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/100);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/150);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/190);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/180);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/170);
sqnorm(URK7(:,end))
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/160);
sqnorm(URK7(:,end))
norm(Uex(:,end)-URK7(:,end))
URK7uf = RK7uf(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/160);
norm(URK7uf-URK7(:, end))
[allNtRK7, allmvRK7, allerRK7] =  error2q2rRK(T2mb2, 160, 20);
URK7 = RK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/160);
sqnorm(URK7(:,end))
[allNtRK7, allmvRK7, allerRK7] =  error2q2rRK(T2mb2, 160, 20);
[allNtRK7, allmvRK7, allerRK7] =  error2q2rRK7(T2mb2, 160, 20);
polyfit(log10(allmvRK7(4:15)), log10(allerRK7(4:15)), 1)
hold on
plot(log10(allmv7), log10(aller7), '-o')
plot(log10(allmv93(2:end)), log10(aller93(2:end)), '-o')
save qubits2HOtest_data
%-- 27/07/2022 9:43 --%
test_harmonic
[U, mniter, matvecs, max_errors, history] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
[U, mniter, matvecs, max_errors, history] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
history
max(history.fUerror)
max(history.f_error)
(factorial(9)*(200/(10*200))^5/(6.7351e-07
))^(1/10)
(factorial(9)*(200/(10*200))^5/(6.7351e-07))^(1/18)
(factorial(9)*(200/(10*200))^9/(6.7351e-07))^(1/18)
200/ans
10/141
[U140, mniter140, matvecs140, max_errors140, history140] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 140, Nt_ts, Ncheb, tol);
norm(U140(:,end))
load errors_harmonic Uex_phase
norm(U(:,end) - Uex_phase)
norm(U(:,end))
norm(Uex_phase)
norm(Uex_phase(:,end))
norm(U(:,end) - Uex_phase(:,end))
norm(U140(:,end) - Uex_phase(:,end))
[U130, mniter130, matvecs130, max_errors130, history130] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 140, Nt_ts, Ncheb, tol, 10, 16, false);
matvecs140
matvecs130
[U130, mniter130, matvecs130, max_errors130, history130] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, Nt_ts, Ncheb, tol, 10, 16, false);
norm(U140(:,end))
norm(U130(:,end))
[U135, mniter135, matvecs135, max_errors135, history135] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 135, Nt_ts, Ncheb, tol, 10, 16, false);
norm(U135(:,end))
load errors_harmonic P
U140E = P\U140;
figure
plot(t, U140E)
plot(t, U140E.*conj(U140E))
plot(t, log10(U140E.*conj(U140E)))
U135E = P\U135;
plot(t, log10(U135E.*conj(U135E)))
[U142, mniter142, matvecs142, max_errors142, history142] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 142, Nt_ts, Ncheb, tol, 10, 16, false);
UE142 = P\UE142;
UE142 = P\U142;
figure
plot(t, log10(U142E.*conj(U142E)))
plot(t, log10(UE142.*conj(UE142)))
clear U142 mniter142 matvecs142 max_errors142 history142 UE142
(factorial(9)*(200/(10*200))^9/(6.7351e-07))^(1/18)
(factorial(9)*(200/(10*200))^9/(6.7351e-06))^(1/18)
200/ans
[U150, mniter150, matvecs150, max_errors150, history150] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, Nt_ts, Ncheb, tol, 10, 16, false);
U150E = P\U150E;
U150E = P\U150;
plot(t, log10(UE150.*conj(UE150)))
plot(t, log10(U150E.*conj(U150E)))
[U160, mniter160, matvecs160, max_errors160, history160] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 160, Nt_ts, Ncheb, tol, 10, 16, false);
U160E = P\U160;
plot(t, log10(U160E.*conj(U160E)))
save Uexact_forced_harmonic Uex
Uex = Uex_phase(:,end)
save Uexact_forced_harmonic Uex
[allNt9, allmv9, aller9, max_ers9] = error_forced_harmonic(9, 9, 140, 15);
allNt9
[allNt9, allmv9, aller9, max_ers9] = error_forced_harmonic(9, 9, 160, 15);
[allNt9, allmv9, aller9, max_ers9] = error_forced_harmonic(9, 9, 170, 10);
polyfit(log10(allmv9(2:5)), log10(aller9(2:5)), 1)
[allNt7, allmv7, aller7, max_ers7] = error_forced_harmonic(7, 7, 140, 10);
[allNt7, allmv7, aller7, max_ers7] = error_forced_harmonic(7, 7, 170, 15);
allNt7(2)
[allNt7, allmv7, aller7, max_ers7] = error_forced_harmonic(7, 7, 200, 12);
polyfit(log10(allmv7(1:8)), log10(aller7(1:8)), 1)
[allNt5, allmv5, aller5, max_ers5] = error_forced_harmonic(5, 5, 200, 15);
[U180M5i, mniter180M5i, matvecs180M5i, max_errors180M5i, history180M5i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 180, 5, 5, eps, 1, 16, false);
[U190M5i, mniter190M5i, matvecs190M5i, max_errors190M5i, history190M5i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 190, 5, 5, eps, 1, 16, false);
norm(U190M5i(:,end))
clear U190M5i mniter190M5i matvecs190M5i max_errors190M5i history190M5i
max_errors190M5i
max_errors180M5i
polyfit(log10(allmv5(1:12)), log10(aller5(1:12)), 1)
[allNt57, allmv57, aller57, max_ers57] = error_forced_harmonic(5, 7, 200, 15);
polyfit(log10(allmv57(1:10)), log10(aller57(1:10)), 1)
figure
plot(log10(allmv9), log10(aller9), '-o')
hold on
plot(log10(allmv7), log10(aller7), '-o')
plot(log10(allmv5), log10(aller5), '-o')
plot(log10(allmv57), log10(aller57), '-o')
[allNt37, allmv37, aller37, max_ers37] = error_forced_harmonic(3, 7, 200, 15);
polyfit(log10(allmv37(1:3)), log10(aller37(1:3)), 1)
polyfit(log10(allmv37(4:end)), log10(aller37(4:end)), 1)
plot(log10(allmv37), log10(aller37), '-o')
[allNt95, allmv95, aller95, max_ers95] = error_forced_harmonic(9, 5, 200, 15);
plot(log10(allmv95(2:end)), log10(aller95(2:end)), '-o')
[allNt93, allmv93, aller93, max_ers93] = error_forced_harmonic(9, 5, 220, 15);
plot(log10(allmv93), log10(aller93), '-o')
log(2)
log(5)
log10(2)
[allNtRK, allmvRK, allerRK] =  error_forced_harmonicRK4(200, 30)
[allNtRK, allmvRK, allerRK] =  error_forced_harmonicRK4(3000, 22)
polyfit(log10(allmvRK(1:11)), log10(allerRK(1:11)), 1)
hold on
plot(log10(allmvRK(1:15)), log10(allerRK(1:15)), '-o')
[allNtRK7, allmvRK7, allerRK7] =  error_forced_harmonicRK7(1e3, 15)
[allNtRK7, allmvRK7, allerRK7] =  error_forced_harmonicRK7(500, 15)
[allNtRK7, allmvRK7, allerRK7] =  error_forced_harmonicRK7(600, 15);
polyfit(log10(allmvRK7(1:7)), log10(allerRK7(1:7)), 1)
hold on
plot(log10(allmvRK7(1:15)), log10(allerRK7(1:15)), '-o')
whos
clear history history130 history135 history140 history150 history160 history180M5i
whos
save test_forced_harmonic
10^0.4
10^0.5
10^0.7
xlabel('log(matvecs)')
ylabel('log(error)')
hold on
clear all
load qubits2HOtest_data
whos
plot(log10(allmvRK7(1:15)), log10(allerRK7(1:15)), '-o')
plot(log10(allmvRK7), log10(allerRK7), '-o')
figure
plot(log10(allmvPWC20), log10(allerPWC20), '-o')
hold on
plot(log10(allmvPWC20Hb), log10(allerPWC20Hb), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
whos
max_errors1
(factorial(9)*(280/(T2mb2*35))^9/(9.2074e-08))^(1/18)
280/ans
max_errorsmax
(factorial(9)*(140/(T2mb2*35))^9/(6.3259e-09))^(1/18)
140/ans
E
allNtRK
Efield = zeros(1, 100);
Efield = zeros(1, 101);
for j=1:101, Efield(:,j) = eig(Hu2modesa); end
Efield = zeros(19, 101);
for j=1:101, Efield(:,j) = eig(Hu2modesa); end
figure
for j=1:101, Efield(:,j) = eig(Hu2modesa - (j-1)*H1umodesa); end
for j=1:101, Efield(:,j) = eig(Hu2modesa - (j-1)*H1u2modesa); end
plot(0:0.01:1, Efield)
clear Efield
Efield1 = zeros(19, 101);
Efield2 = zeros(19, 101);
for j=1:101, Efield(:,j) = eig(Hu2modesa - (j-1)*0.01*H1u2modesa); end
figure
plot(0:0.01:1, Efield1)
for j=1:101, Efield1(:,j) = eig(Hu2modesa - (j-1)*0.01*H1u2modesa); end
clear Efield
for j=1:101, Efield2(:,j) = eig(Hu2modesa - (j-1)*0.01*H1u2modesa); end
plot(0:0.01:1, Efield1)
figure
plot(0:0.01:1, log10(abs(Efield1-Efield(:,1))))
plot(0:0.01:1, log10(abs(Efield1-Efield1(:,1))))
plot(0:0.01:1, log10(abs((Efield1-Efield1(:,1))./Efield1(:,1))))
E
iEwithout0 = [1:6, 10:19];
plot(0:0.01:1, log10(abs((Efield1(iEwithout0,:)-Efield1(iEwithout0,1))./Efield1(iEwithout0,1))))
figure
plot(0:0.01:1, log10(abs((Efield1(7:9,:)-Efield1(7:9,1)))))
figure
plot(0:0.01:1, log10(abs((Efield1(7:9,:)-Efield1(7:9,1)))))
plot(0:0.01:1, (((Efield1(7:9,:)-Efield1(7:9,1)))))
figure
plot(0:0.01:1, (((Efield2(7:9,:)-Efield2(7:9,1)))))
for j=1:101, Efield2(:,j) = eig(Hu2modesa - (j-1)*0.01*H2u2modesa); end
figure
plot(0:0.01:1, (((Efield2(7:9,:)-Efield2(7:9,1)))))
figure
plot(0:0.01:1, log10(abs((Efield2(iEwithout0,:)-Efield2(iEwithout0,1))./Efield2(iEwithout0,1))))
xlabel('$\frac{f_1(t)}{g_2}$', 'interpreter', 'latex')
max(abs(allfieldmax), [], 2)
ylabel('$E_n$', 'interpreter', 'latex')
figure
plot(0:0.01:1, Efield1)
clf
plot(0:0.01:1, Efield1(iEwithout0, :))
figure
plot(0:0.01:1, Efield1(7:9, :))
xlabel('$\frac{f_1(t)}{g_2}$', 'interpreter', 'latex')
ylabel('$E_n$', 'interpreter', 'latex')
xlabel('$\frac{f_1(t)}{g_2}$', 'interpreter', 'latex')
ylabel('$E_n$', 'interpreter', 'latex')
xlabel('$\frac{f_1(t)}{g_2}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left\frac{|E_n(f_1)-E_n(0)\right|}{\left|E_n(0)\right|}\right)$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\frac{|E_n(f_1)-E_n(0)|}{|E_n(0)|}\right)$', 'interpreter', 'latex')
xlabel('$\frac{f_1}{g_2}$', 'interpreter', 'latex')
save qubits2HOtest_data
figure
plot(0:pi/T2mb2:50, fieldw2mb2)
[H2modesa60, H1_2modesa60, H2_2modesa60, Hu2modesa60, H1u2modesa60, H2u2modesa60, Hoperations2modesa60, fcouplingOp2modesa60] = generateH2mode_gen(0.2/g2GHz, 1/sqrt(2), 1, 60);
[P60, D60] = eig(full(Hu2modesa60));
E60 = diag(D60);
E60
whos
eig(Hu2modesa60 + 1.07*H1u2modesa60 + 1.07*H2u2modesa60)
eig(Hu2modesa60 - 1.07*H1u2modesa60 - 1.07*H2u2modesa60)
60*0.03
[UDw60, mnitercDw60, matvecsDw60, max_errorsDw60, historyDw60] = SemiGlobalHparams(Hoperations2modesa60.psi, Hoperations2modesa60.diff_psi, 0, allfield, [], [-15 121], u02modes, 0:T2mb2/280:T2mb2, Nt, 9, 9, 1e-6, 10, 16, [], false);
max_errorsDw60
norm(UDw60(:,end) - U1(:,end))
figure
UDw60(:,end).*conj(UDw60(:,end))
sqnorm(UDw60(:,end) - U1(:,end))/4
max_errorsDw60
(factorial(9)*(140/(T2mb2*120))^9/(9.7253e-08))^(1/18)
mniterDw60
mnitercDw60
mniterc1
Gop2mb2_60 = @(u, t, v) -1i*(Hu2modesa60*v - [H1u2modesa60*v, H2u2modesa60*v]*get_fieldt(t))
Gdiff_op2mb2_60 = @(u1, t1, u2, t2) 1i*((get_fieldt1(t1) - get_fieldt1(t2)).*(H1u2modesa60*u1) + (get_fieldt2(t1) - get_fieldt2(t2)).*(H2u2modesa60*u1))
280/1.2860
[UDw60_220, mnitercDw60_220, matvecsDw60_220, max_errorsDw60_220, historyDw60_220] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/220:T2mb2, 220, 9, 9, 1e-6, 10, 16, false);
mnitercDw60_220
sqnorm(UDw60_220(:,end))
max_errorsDw60_220
(factorial(9)*(280/(T2mb2*120))^9/(9.7253e-08))^(1/18)
280/ans
(factorial(9)*(220/(T2mb2*120))^9/(8.0907e-07))^(1/18)
220/ans
clear UDw60_220 mnitercDw60_220 matvecsDw60_220 max_errorsDw60_220 historyDw60_220
[UDw60_155, mnitercDw60_155, matvecsDw60_155, max_errorsDw60_155, historyDw60_155] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/155:T2mb2, 155, 9, 9, 1e-6, 10, 16, false);
sqnorm(UDw60_155(:,end))
max_errorsDw60_155
(factorial(9)*(155/(T2mb2*120))^9/(1.6146e-05))^(1/18)
matvecsDw60_155
matvecsDw60
mnitercDw60_155
[UDw60_150, mnitercDw60_150, matvecsDw60_150, max_errorsDw60_150, historyDw60_150] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/150:T2mb2, 150, 9, 9, 1e-6, 10, 16, false);
sqnorm(UDw60_150(:,end))
max_errorsDw60_150
[UDw60_145, mnitercDw60_145, matvecsDw60_145, max_errorsDw60_145, historyDw60_145] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/150:T2mb2, 145, 9, 9, 1e-6, 10, 16, false);
sqnorm(UDw60_150(:,end))
sqnorm(UDw60_145(:,end))
[UDw60_140, mnitercDw60_140, matvecsDw60_140, max_errorsDw60_140, historyDw60_140] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/140:T2mb2, 140, 9, 9, 1e-6, 10, 16, false);
sqnorm(UDw60_140(:,end))
max_errorsDw60_140
(factorial(9)*(140/(T2mb2*120))^9/(max_errorsDw60_140.f))^(1/18)
[UDw60_135, mnitercDw60_135, matvecsDw60_135, max_errorsDw60_135, historyDw60_135] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/135:T2mb2, 135, 9, 9, 1e-6, 10, 16, false);
sqnorm(UDw60_135(:,end))
mnitercDw60_140
mnitercDw60_145
[UDw60_145, mnitercDw60_145, matvecsDw60_145, max_errorsDw60_145, historyDw60_145] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/145:T2mb2, 145, 9, 9, 1e-6, 10, 16, false);
plot(0:T2mb2/140:T2mb2, log10(UDw60_145E.*conj(UDw60_145E)))
UDw60_145E = P60\UDw60_145;
plot(0:T2mb2/140:T2mb2, log10(UDw60_145E.*conj(UDw60_145E)))
plot(0:T2mb2/140:T2mb2, log10(UDw60_140E.*conj(UDw60_140E)))
plot(0:T2mb2/145:T2mb2, log10(UDw60_145.*conj(UDw60_145)))
plot(0:T2mb2/145:T2mb2, log10(UDw60_145E.*conj(UDw60_145E)))
UDw60_140E = P60\UDw60_140;
UDw60_135E = P60\UDw60_135;
plot(0:T2mb2/140:T2mb2, log10(UDw60_140E.*conj(UDw60_140E)))
figure
UDw60E = P60\UDw60;
plot(0:T2mb2/280:T2mb2, log10(UDw60_280E.*conj(UDw60_280E)))
plot(0:T2mb2/280:T2mb2, log10(UDw60E.*conj(UDw60E)))
figure
plot(0:T2mb2/135:T2mb2, log10(UDw60_135E.*conj(UDw60_135E)))
matvecsDw145
matvecsDw60_145
[UDw60_145_2, mnitercDw60_145_2, matvecsDw60_145_2, max_errorsDw60_145_2, historyDw60_145_2] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/145:T2mb2, 145, 9, 9, 1e-3, 10, 16, false);
mnitercDw60_145_2
max_errorsDw60_145_2
sqnorm(UDw60_145_2(:,end))
norm(UDw60_145_2(:,end) - UDw60_145_2(:,end))
norm(UDw60_145_2(:,end) - UDw60(:,end))
[UDw60_145_3, mnitercDw60_145_3, matvecsDw60_145_3, max_errorsDw60_145_3, historyDw60_145_3] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/145:T2mb2, 145, 9, 9, 1e-1, 10, 16, false);
max_errorsDw60_145_3
mnitercDw60_145_3
matvecsDw60_145_3
norm(UDw60_145_3(:,end) - UDw60(:,end))
clear UDw60_145_2 mnitercDw60_145_2 matvecsDw60_145_2 max_errorsDw60_145_2 historyDw60_145_2
clear UDw60_155 mnitercDw60_155 matvecsDw60_155 max_errorsDw60_155 historyDw60_155
clear UDw60_150 mnitercDw60_150 matvecsDw60_150 max_errorsDw60_150 historyDw60_150
[UDw60_145_3, mnitercDw60_145_3, matvecsDw60_145_3, max_errorsDw60_145_3, historyDw60_145_3] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/145:T2mb2, 145, 9, 9, 1e-1, 10, 16, false);
sqnorm(Unew)
sqnorm(Uguess)
(factorial(5)*(140/(T2mb2*120))^5/(1e-4))^(1/18)
(factorial(5)*(140/(T2mb2*120))^5/(1e-4))^(1/10)
[UDw60_150M5, mnitercDw60_150M5, matvecsDw60_150M5, max_errorsDw60_150M5, historyDw60_150M5] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/150:T2mb2, 150, 5, 5, 1e-1, 10, 16, false);
sqnorm(UDw60_150M5(:,end))
max_errorsDw60_150M5
(factorial(5)*(140/(T2mb2*120))^5/(4.8e-2))^(1/10)
150/ans
whos
save qubits2HOtest_data
%-- 29/07/2022 0:22 --%
load qubits2HOtest_data
(factorial(9)*(140/(T2mb2*120))^9/(max_errorsDw60_140.f))^(1/18)
140/ans
Nt
[UDw60, mnitercDw60, matvecsDw60, max_errorsDw60, historyDw60] = SemiGlobalHparams(Hoperations2modesa60.psi, Hoperations2modesa60.diff_psi, 0, allfield, [], [-15 121], u02modes, 0:T2mb2/280:T2mb2, Nt, 9, 9, 1e-6, 10, 16, [], false);
[UDw60_155, mnitercDw60_155, matvecsDw60_155, max_errorsDw60_155, historyDw60_155] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/155:T2mb2, 155, 9, 9, 1e-6, 10, 16, false);
vecnorm(Ulast)
whos
[U43, mniterc43, matvecs43, max_errors43, history43] = SemiGlobal(Gop2mb2, Gdiff_op2mb2, 0, [], [-35*1i 15*1i], u02modes, 0:T2mb2/43:T2mb2, 43, 9, 9, 1e-6, 10, 16, false);
vecnorm(Ulast)
figure
plot(0:T2mb2/135:T2mb2, log10(UDw60_135E.*conj(UDw60_135E)))
figure
plot(0:T2mb2/280:T2mb2, log10(UDw60E.*conj(UDw60E)))
figure
plot(0:T2mb2/140:T2mb2, log10(UDw60_140E.*conj(UDw60_140E)))
plot((0:T2mb2/280:T2mb2)/(2*pi), log10(UDw60E.*conj(UDw60E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
figure
plot((0:T2mb2/140:T2mb2)/(2*pi), log10(UDw60_140E.*conj(UDw60_140E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
title('$\Delta\omega=60,\Delta t = T/140$', 'interpreter', 'latex')
plot((0:T2mb2/135:T2mb2)/(2*pi), log10(UDw60_135E.*conj(UDw60_135E)))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_n|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
title('$\Delta\omega=60,\Delta t = T/135$', 'interpreter', 'latex')
whos
b8 = 77/1440;
weightsRK7 = [77/1440 - b8; 0; 0; 32/105; 1771561/6289920; 243/2560; 16807/74880; b8; 11/270];
alphaT7 = [1/6                       0       0           0               0                   0               0                   0;
0                         1/3     0           0               0                   0               0                   0;
1/8                       0       3/8         0               0                   0               0                   0;
148/1331                  0       150/1331    -56/1331        0                   0               0                   0;
-404/243                  0       -170/27     4024/1701       10648/1701          0               0                   0;
2466/2401                 0       1242/343    -19176/16807    -51909/16807        1053/2401       0                   0;
1/(576*b8)                0       0           1/(105*b8)      -1331/(279552*b8)   -9/(1024*b8)    343/(149760*b8)     0;
-71/32 - 270*b8/11        0       -195/22     32/7            29403/3584          -729/512        1029/1408           270*b8/11].';
cRK7 = [1/6, 1/3, 1/2, 2/11, 2/3, 6/7, 0, 1];
URK7gen = RKgeneral(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/160, cRK7, weightsRK7, alphaT7);
max(abs(URK7gen(:,end)- URK7(:,end))
max(abs(URK7gen(:,end)- URK7(:,end)))
URK7gen_uf = RKgeneral_uf(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)),  [0 T2mb2], u02modes, T2mb2/160, cRK7, weightsRK7, alphaT7);
max(abs(URK7gen_uf- URK7(:,end)))
clear URK7gen URK7gen alphaT7 cRK7 b8 weightsRK7
[allNtRK7a, allmvRK7a, allerRK7a] =  error_decayRK7(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)), [0 T2mb2], u02modes, Uex(:,end), 160, 20);
allerRK7 - allerRK7a
clear allNtRK7a allmvRK7a allerRK7a
whos
[UDw60_1e3, mnitercDw60_1e3, matvecsDw60_1e3, max_errorsDw60_1e3] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/1e3:T2mb2, 1e3, 9, 13, eps, 10, 16, false);
max_errorsDw60_1e3
[UDw60_1200, mnitercDw60_1200, matvecsDw60_1200, max_errorsDw60_1200] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/1200:T2mb2, 1200, 9, 13, eps, 10, 16, false);
norm(UDw60_1200(:,end) - UDw60_1e3(:,end))
norm(UDw60_1200(:,end) - UDw60(:,end))
clear UDw60_1e3 mnitercDw60_1e3 matvecsDw60_1e3 max_errorsDw60_1e3
[allNt9Dw60, allmv9Dw60, aller9Dw60, max_ers9Dw60] = error_decaySG(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, UDw60_1200(:,end),...
[0 T2mb2], 9, 9, 180, 15, 1);
[allNt9Dw60, allmv9Dw60, aller9Dw60, max_ers9Dw60] = error_decaySG(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, UDw60_1200(:,end), [0 T2mb2], 9, 9, 160, 15, 1);
[allNt9Dw60, allmv9Dw60, aller9Dw60, max_ers9Dw60] = error_decaySG(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, UDw60_1200(:,end), [0 T2mb2], 9, 9, 150, 15, 1);
[UDw60_145i, mnitercDw60_145i, matvecsDw60_145i, max_errorsDw60_145i] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/145:T2mb2, 145, 9, 9, eps, 1, 16, false);
sqnorm(UDw60_145i(:,end))
[UDw60_140i, mnitercDw60_140i, matvecsDw60_140i, max_errorsDw60_140i] = SemiGlobal(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, 0:T2mb2/140:T2mb2, 140, 9, 9, eps, 1, 16, false);
sqnorm(UDw60_140i(:,end))
[allNt9Dw60, allmv9Dw60, aller9Dw60, max_ers9Dw60] = error_decaySG(Gop2mb2_60, Gdiff_op2mb2_60, 0, [], [-121*1i 15*1i], u02modes, UDw60_1200(:,end), [0 T2mb2], 9, 9, 150, 15, 1);
polyfit(log10(allmv9Dw60(3:12)), log10(aller9Dw60(3:12)), 1)
[allNt, allmv, aller] =  error_decayRK4(@(t, u) -1i*(Hu2modesa60*u - [H1u2modesa60*u, H2u2modesa60*u]*get_fieldt(t)), [0 T2mb2], u02modes, UDw60_1200(:, end), 300, 15);
[allNt, allmv, aller] =  error_decayRK4(@(t, u) -1i*(Hu2modesa60*u - [H1u2modesa60*u, H2u2modesa60*u]*get_fieldt(t)), [0 T2mb2], u02modes, UDw60_1200(:, end), 300, 30);
hold on
plot(log10(allmv9Dw60), log10(aller9Dw60), '-o')
save qubits2HOtest_data
%-- 04/08/2022 10:20 --%
load qubits2HOtest_data
whos
E60
u02modesa
u02modes
P60\u02modes
(P60\u02modes).*conj((P60\u02modes))
figure
plot((0:T2mb2/280:T2mb2)/(2*pi), (UDw60E(19,:).*conj(UDw60E(19,:))) - (UDw60E(19,1).*conj(UDw60E(19,1))))
UDw60_1200E = P60\UDw60_1200;
plot((0:T2mb2/280:T2mb2)/(2*pi), (UDw60_1200E(19,:).*conj(UDw60_1200E(19,:))) - (UDw60_1200E(19,1).*conj(UDw60_1200E(19,1))))
plot((0:T2mb2/1200:T2mb2)/(2*pi), (UDw60_1200E(19,:).*conj(UDw60_1200E(19,:))) - (UDw60_1200E(19,1).*conj(UDw60_1200E(19,1))))
figure
plot((0:T2mb2/280:T2mb2)/(2*pi), log10(UDw60_1200E(19,:).*conj(UDw60_1200E(19,:))))
plot((0:T2mb2/1200:T2mb2)/(2*pi), log10(UDw60_1200E(19,:).*conj(UDw60_1200E(19,:))))
plot((0:T2mb2/1200:T2mb2)/(2*pi), (UDw60_1200E(19,:).*conj(UDw60_1200E(19,:))))
xlabel('$\frac{g_2 t}{2\pi}$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(\left|\left<\varphi_{19}|\psi(T)\right>\right|^2\right)$', 'interpreter', 'latex')
ylabel('$\left|\left<\varphi_{19}|\psi(T)\right>\right|^2$', 'interpreter', 'latex')
(UDw60_1200E(19,end).*conj(UDw60_1200E(19,end))) - (UDw60_1200E(19,1).*conj(UDw60_1200E(19,1)))
title('$\Delta\omega=60,\Delta t = T/1200$', 'interpreter', 'latex')
E60(end)
UDw60_1200E(:,1).*conj(UDw60_1200E(:,1))
[Pd60, Dd60] = eig(Hu2modesa(10:19, 10:19));
[Pd60, Dd60] = eig(full(Hu2modesa(10:19, 10:19)));
Ed60 = diag(Dd60)
[Pd60, Dd60] = eig(full(Hu2modesa60(10:19, 10:19)));
Ed60 = diag(Dd60)
Pd60\u02modes(10:19)
[Ed60, abs(Pd60\u02modes(10:19))^2]
[Ed60, abs(Pd60\u02modes(10:19)).^2]
[Ed60, abs(Pd60\UDw60_1200(10:19, end)).^2]
%-- 09/08/2022 10:17 --%
load qubits2HOtest_data
figure
whos
plot(log10(allNtPWC20), log10(allerPWC20), '-o')
hold on
plot(log10(allNtPWC20Hb), log10(allerPWC20Hb), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')