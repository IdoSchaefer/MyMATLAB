test_harmonicSG1
Uex_forced_harmonic
load Uexact_forced_harmonic
max(abs(Uex - Uex_fHO))
test_harmonicSG1
[U0, mniter0, matvecs0, max_errors0, history0] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
matvecs0
mniter0
max_errors0
[U, mniter, matvecs, max_errors, history] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
est_errors
doc norm
error_anal_t = vecnorm(U - Uex_phase)./vecnorm(Uex_phase);
figure
plot(0:10/150:10, error_anal_t)
size(error_anal_t)
size(t)
plot(0:0.1:10, error_anal_t)
hold on
plot(2/3:2/3:10, history.fUerror)
size(history.fUerror)
size(2/3:2/3:10)
size(1/15:1/15:10)
plot(1/15:1/15:10, history.fUerror)
plot(1/15:1/15:10, cumsum(history.fUerror))
[U1, mniter1, matvecs1, est_errors1, history1] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
error_anal1 = norm(U1(:, end) - Uex_fHO)/norm(Uex_fHO)
mniter
matvecs
est_errors1
figure
error_anal_t1 = vecnorm(U1(:, end) - Uex_fHO)./vecnorm(Uex_fHO);
figure
plot(0:0.1:10, error_anal_t1)
error_anal_t1 = vecnorm(U1(:, end) - Uex_phase)./vecnorm(Uex_phase);
plot(0:0.1:10, error_anal_t1)
error_anal_t1 = vecnorm(U1 - Uex_phase)./vecnorm(Uex_phase);
plot(0:0.1:10, error_anal_t1)
hold on
plot(1/20:1/20:10, cumsum(history1.fUerror))
[U1o, mniter1o, matvecs1o, est_errors1o, history1o] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
[mniter1o, matvecs1o]
error_anal_t1o = vecnorm(U1o - Uex_phase)./vecnorm(Uex_phase);
figure
plot(0:0.1:10, error_anal_t1)
hold on
plot(0:0.1:10, error_anal_t1o)
lambda = conv_ratios_1st(p_max)
lambda = conv_ratios_1st(10)
(factorial(9)*(140/(T2mb2*35))^9/(6.326e-9))^(1/18)
[U1o, mniter1o, matvecs1o, est_errors1o, history1o] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
[U1, mniter1, matvecs1, est_errors1, history1] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
[Um6, mniterm6, matvecsm6, est_errorsm6, historym6] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 6, 7, tol);
tcheb6 = -cos(((1:6) - 1)*pi/(6-1));
tcheb6
acos(0)
[Um6, mniterm6, matvecsm6, est_errorsm6, historym6] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 6, 7, tol);
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
max(max(abs(U1-U1test)))
est_errors1
est_errors1test
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
max(max(abs(U1-U1test)))
est_errors1test
est_errors1
est_errors1test.texp_error_cheap/est_errors1.texp_error_cheap
est_errors1test.texp_error_exact/est_errors1.texp_error_exact
est_errors1test.texp_exact/est_errors1.texp_exact
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
est_errors1test.texp_exact/est_errors1.texp_exact
est_errors1test.texp_error_exact/est_errors1.texp_error_exact
est_errors1test.texp_error_cheap/est_errors1.texp_error_cheap
est_errors1test.texp_cheap/est_errors1.texp_cheap
est_errors1test.texp_cheap/est_errors1.texp_error_cheap
0.05*195
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
195^5*f_scalar_error
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
195^5*f_scalar_error/factorial(5)
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
cos(Nfm*fm_cheb_samplingp)
cos((Nfm+1)*fm_cheb_samplingp)
cos((Nfm)*acos(fm_cheb_samplingp))
fm_sp
abs(dvd_fz4)
f_fun(fm_sp, Tts)
dvd_fz4
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
abs(dvd_fz4)
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
est_errors1
est_errors1test
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
est_errors1test
est_errors1
4*3.5901e-06
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
est_errors1
est_errors1test
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
est_errors1test
est_errors1
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 5, 7, tol);
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 9, 9, tol);
norm(U1test(:,end))
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 120, 9, 9, tol);
norm(U1test(:,end))
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, tol);
norm(U1test(:,end))
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 140, 9, 9, tol);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, tol);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
t
[~, ~, ~, E, P, H] = gsV(V + x, xdomain, Nx);
[~, ~, ~, E, P, H] = gsV(V + x, [-L/2, L/2], Nx);
E
[~, ~, ~, E, P, H] = gsV(V - x, [-L/2, L/2], Nx);
[~, ~, ~, E, P, H] = gsV(V, [-L/2, L/2], Nx);
E
[U1test, mniter1test, matvecs1test, est_errors1test, history1test] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, tol);
fztest
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
[~, ~, ~, E, P, H] = gsV(V, [-L/2, L/2], Nx);
U1testE = P\U1test;
figure
plot(0:1/13:10, U1testE.*conj(U1testE))
size(U1test)
plot(0:0.1:10, U1testE.*conj(U1testE))
plot(0:0.1:10, log10(U1testE.*conj(U1testE)))
figure
plot(t, cos(t))
[U130, mniter130, matvecs130, est_errors130, history130] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, tol);
[U120, mniter120, matvecs120, est_errors120, history120] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 120, 9, 9, tol);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts>1;
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts>
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
U130 = P\U130;
U120 = P\U120;
[U140, mniter140, matvecs140, est_errors140, history140] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 140, 9, 9, tol);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
[U120, mniter120, matvecs120, est_errors120, history120] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 120, 9, 9, tol);
[U130, mniter130, matvecs130, est_errors130, history130] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, tol);
U130E = P\U130;
U120E = P\U120;
U140E = P\U140;
figure
plot(0:0.1:10, log10(U140E.*conj(U140E)))
norm(U1test(:,end))
norm(U130(:,end))
norm(U130(:,end))^2
39/130
20/130
20^(1/130)
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
figure
U150E = P\U150;
plot(0:0.1:10, log10(U150E.*conj(U150E)))
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
max(abs(ev_domain))
f_scalar_error/max(abs(fztest))
f_scalar_error/min(abs(fztest))
4.6603e-06
4.6603e-06*(195*1/15)^Nt_ts/factorialNt_ts
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
f_scalar_error_rel = max(abs(chebc2result(Ccheb_f_ts(:, Nt_ts - 1), ev_domain, ztest) - fztest)./abs(fztest));
max(abs(chebc2result(Ccheb_f_ts(:, Nt_ts - 1), ev_domain, ztest) - fztest)./abs(fztest))
195/15
13^Nt_ts/factorialNt_ts
tol
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
figure
U150Ea = P\U150a;
plot(0:0.1:10, log10(U150Ea.*conj(U150Ea)))
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
f_scalar_error
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
f_scalar_error
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8)
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
195^9/factorial(9)
f_fun(-195*1i, 1/15, 9, 1e-5/150)
f_fun(-195*1i, 1/15, 9, 1e-8/150)
f_fun(-195*1i, 1/15, 9, 1e-5/150)
f_fun(-195*1i, 1/15, 9, 1e-5/150)-f_fun(-195*1i, 1/15, 9, 1e-8/150)
f_fun(-1i, 1/15, 9, 1e-5/150)-f_fun(-1i, 1/15, 9, 1e-8/150)
mniter150
mniter150a
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
1e-8/150
[U150b, mniter150b, matvecs150b, est_errors150b, history150b] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
figure
U150Eb = P\U150b;
plot(0:0.1:10, log10(U150Eb.*conj(U150Eb)))
figure
plot(0:0.1:10, log10(U150E.*conj(U150E)))
[U150b, mniter150b, matvecs150b, est_errors150b, history150b] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-10);
U150Eb = P\U150b;
figure
plot(0:0.1:10, log10(U150Eb.*conj(U150Eb)))
[U150b, mniter150b, matvecs150b, est_errors150b, history150b] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-10);
[U150a, mniter150a, matvecs150a, est_errors150a, history150a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-8);
mniter150a
mniter150b
history150.niter
[history150.niter; history150a.niter; history150b.niter; 0:1/15:10]
[history150.niter; history150a.niter; history150b.niter; 1/15:1/15:10]
[U150h, mniter150h, matvecs150h, est_errors150h, history150h] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + 0.5*x*cos(t), v), @(u1, t1, u2, t2) -1i*x*0.5*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-5);
U150Eh = P\U150h;
figure
plot(0:0.1:10, log10(U150Eh.*conj(U150Eh)))
10^-7.5
1e-5/150
figure
plot(0:0.1:10, log10(U140.*conj(U140)))
plot(0:0.1:10, log10(U140E.*conj(U140E)))
figure
plot(0:0.1:10, log10(U130E.*conj(U130E)))
[U130a, mniter130a, matvecs130a, est_errors130a, history130a] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, 1e-8);
U130aE = P\U130a;
figure
plot(0:0.1:10, log10(U130Ea.*conj(U130Ea)))
plot(0:0.1:10, log10(U130aE.*conj(U130aE)))
[U130b, mniter130b, matvecs130b, est_errors130b, history130b] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, 1e-10);
U130bE = P\U130b;
figure
plot(0:0.1:10, log10(U130bE.*conj(U130bE)))
[U130o, mniter130o, matvecs130o, est_errors130o, history130o] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, 1e-8);
mniter130o
mniter130
U130oE = P\U130o;
[U130o, mniter130o, matvecs130o, est_errors130o, history130o] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, 1e-5);
mniter130
mniter130o
norm(U130o(:,end))^2
norm(U130(:,end))^2
150*0.9
1e12^(1/135)
1e9^(1/75)
8.5/5
1e8.5^(1/75)
(10^8.5)^(1/75)
log10(e)
e
log10(exp(1))
1.7/ans
exp(3.9/150)
exp(3.9/15)
[U130b, mniter130b, matvecs130b, est_errors130b, history130b] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, 1e-10);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
[U150h, mniter150h, matvecs150h, est_errors150h, history150h] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + 0.5*x*cos(t), v), @(u1, t1, u2, t2) -1i*x*0.5*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-5);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
xdomain
L/2
L/2*(1/15)^2/20
[U140, mniter140, matvecs140, est_errors140, history140] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 140, 9, 9, tol);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
max(abs(chebc2result(Ccheb_f_next(:, Nt_ts - 1), ev_domain, ztest) - f_fun(ztest, 2*Tts)))
size(Ccheb_f_next)
max(abs(chebc2result(Ccheb_f_next(:, Nt_ts - 1), ev_domain, ztest) - f_fun(ztest, 2*Tts))./abs(f_fun(ztest, 2*Tts)))
abs(f_fun(ztest, 2*Tts)
abs(f_fun(ztest, 2*Tts))
abs(f_fun(ztest, Tts))
figure
plot(1/13:1/13:10, history130.texp_error)
history130
plot(1/13:1/13:10, history130.texp_error.exact)
plot(1/13:1/13:10, history130.reldif)
[U130, mniter130, matvecs130, est_errors130, history130] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 130, 9, 9, tol);
all_gamma = conv_ratios_texp(9, 4)
all_eta = conv_ratios_fm(18, 4)
L/2*(1/13)^2/20
L/2*(1/13)^2.*all_gamma(:,end)
5.821174e+00*3.06e-3
7.959365e-03*2.4361e-03
[U150h, mniter150h, matvecs150h, est_errors150h, history150h] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + 0.5*x*cos(t), v), @(u1, t1, u2, t2) -1i*x*0.5*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, 1e-5);
[U150, mniter150, matvecs150, est_errors150, history150] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 150, 9, 9, tol);
[U140, mniter140, matvecs140, est_errors140, history140] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 140, 9, 9, tol);
[U140kr, mniter140kr, matvecs140kr, est_errors140kr, history140kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 140, 9, 9, tol);
norm(U140kr(:, end))
1-norm(U140kr(:, end))
est_errors140kr
whos
norm(U140kr(:, end)-Uex)
est_errors140kr
U150krE = P\U150kr;
U140krE = P\U140kr;
figure
plot(0:0.1:10, log10(U140krE.*conj(U140krE)))
[U130kr, mniter130kr, matvecs130kr, est_errors130kr, history130kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 130, 9, 9, tol);
est_errors130kr
norm(U130kr(:, end)-Uex)
[U120kr, mniter120kr, matvecs120kr, est_errors120kr, history120kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 120, 9, 9, tol);
est_errors120kr
mniter120kr
mniter130kr
norm(U120kr(:, end)-Uex)
[U100kr, mniter100kr, matvecs100kr, est_errors100kr, history100kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 100, 9, 9, tol);
est_errors100kr
norm(U100kr(:, end)-Uex)
U100krE = P\U100kr;
plot(0:0.1:10, log10(U100krE.*conj(U100krE)))
[U80kr, mniter80kr, matvecs80kr, est_errors80kr, history80kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 80, 9, 9, tol);
est_errors80kr
norm(U80kr(:, end)-Uex)
mniter80kr
U80krE = P\U80kr;
figure
plot(0:0.1:10, log10(U80krE.*conj(U80krE)))
plot(0:0.1:10, log10(sqrt(U80krE.*conj(U80krE))))
[U90kr, mniter90kr, matvecs90kr, est_errors90kr, history90kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 90, 9, 9, tol);
est_errors90kr
U90krE = P\U90kr;
figure
plot(0:0.1:10, log10(sqrt(U90krE.*conj(U90krE))))
mniter90kr
[U70kr, mniter70kr, matvecs70kr, est_errors70kr, history70kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 70, 9, 9, tol);
norm(U70kr(:, end)-Uex)
U70krE = P\U70kr;
figure
plot(0:0.1:10, log10(sqrt(U70krE.*conj(U70krE))))
[U70kr, mniter70kr, matvecs70kr, est_errors70kr, history70kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 70, 9, 9, tol);
f_scalar_error*max(abs(eigval))^Nt_ts/factorialNt_ts
[U70kr, mniter70kr, matvecs70kr, est_errors70kr, history70kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 70, 9, 9, tol);
max(abs(eigval))
[U100kr, mniter100kr, matvecs100kr, est_errors100kr, history100kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 100, 9, 9, tol);
[U120kr, mniter120kr, matvecs120kr, est_errors120kr, history120kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 120, 9, 9, tol);
[U130kr, mniter130kr, matvecs130kr, est_errors130kr, history130kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 130, 9, 9, tol);
figure
plot(0:0.1:10, log10(sqrt(U120krE.*conj(U120krE))))
U120krE = P\U120kr;
plot(0:0.1:10, log10(sqrt(U120krE.*conj(U120krE))))
figure
plot(0:0.1:10, log10(sqrt(U100krE.*conj(U100krE))))
[U70kr, mniter70kr, matvecs70kr, est_errors70kr, history70kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 70, 9, 9, tol);
f_scalar_error*max(abs(eigval))^Nt_ts/factorialNt_ts
[U100kr, mniter100kr, matvecs100kr, est_errors100kr, history100kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 100, 9, 9, tol);
f_scalar_error*max(abs(eigval))^Nt_ts/factorialNt_ts
[U90kr, mniter90kr, matvecs90kr, est_errors90kr, history90kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 90, 9, 9, tol);
f_scalar_error*max(abs(eigval))^Nt_ts/factorialNt_ts
[U200_97, mniter200_97, matvecs200_97, est_errors200_97, history200_97] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 9, 7, tol);
norm(U200_97(:, end)-Uex)
[U200_95, mniter200_95, matvecs200_95, est_errors200_95, history200_95] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 9, 5, tol);
[U200_96, mniter200_96, matvecs200_96, est_errors200_96, history200_96] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 200, 9, 6, tol);
est_errors200_96
norm(U200_96(:, end)-Uex)
U200_96E = P\U200_96;
figure
plot(0:0.1:10, log10(sqrt(U200_96E.*conj(U200_96E))))
[U200_95kr, mniter200_95kr, matvecs200_95kr, est_errors200_95kr, history200_95kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 200, 9, 5, tol);
norm(U200_95kr(:, end)-Uex)
est_errors200_95kr
[U200_94kr, mniter200_94kr, matvecs200_94kr, est_errors200_94kr, history200_94kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 200, 9, 4, tol);
U200_94krE = P\U200_94;
U200_94krE = P\U200_94kr;
figure
plot(0:0.1:10, log10(sqrt(U200_94krE.*conj(U200_94krE))))
est_errors200_94kr
norm(U200_94kr(:, end)-Uex)
U200_95krE = P\U200_95kr;
figure
plot(0:0.1:10, log10(sqrt(U200_95krE.*conj(U200_95krE))))
[U200_93kr, mniter200_93kr, matvecs200_93kr, est_errors200_93kr, history200_93kr] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, 200, 9, 3, tol);
est_errors200_93kr
U200_93krE = P\U200_93kr;
figure
plot(0:0.1:10, log10(sqrt(U200_93krE.*conj(U200_93krE))))
mniter200_94kr
mniter200_93kr
log10(1e-5/200)
%-- 09/11/2023 13:05 --%
load Uex_article
whos
load coulomb_optV240
whos
Gop = @(u,t,v) -1i*Hpsi(K240, Vabs240 - xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), v)
Gdiff_op=@(u1, t1, u2, t2) 1i*xabs240*0.1*(sech((t1-500)/(170)).^2.*cos(0.06*(t1-500)) - sech((t2-500)/(170))^2*cos(0.06*(t2-500))).*u1
[u, ~, matvecs, all_est_er] = SemiGlobal(@(u, t, v) -1i*Hpsi(K240, Vabs240 - xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), v),...
@(u1, t1, u2, t2) 1i*xabs240*0.1*(sech((t1-500)/(170)).^2.*cos(0.06*(t1-500)) - sech((t2-500)/(170))^2*cos(0.06*(t2-500))).*u1,...
0, [], [], fi0240, [0 T], Nt, Nt_ts, Nkr, eps, Niter, 20, false);
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, [0 1e3], 5e3, 9, 7, 1e-5, 10, 16, false);
est_errors
mniter
whos
norm(U(:, end)-Uex(:,2))
norm(U(:, end)-Uex(:,2))/norm(Uex(:,2))
[U95, mniter95, matvecs95, est_errors95, history95] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, [0 1e3], 5e3, 9, 5, 1e-5, 10, 16, false);
est_errors95
norm(U95(:, end)-Uex(:,2))/norm(Uex(:,2))
mniter
mniter59
mniter95
whos
[~, ~, ~, Etest, P240] = gsV(Vabs240, xdomain240, 768);
max(abs(E-Etest))
max(abs(E240-Etest))
Etest(1:10)
E240(1:10)
E240-Etest
E240
whos
U95E = P240\U95;
figure
[U95, mniter95, matvecs95, est_errors95, history95] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 5e3, 9, 5, 1e-5, 10, 16, false);
U95E = P240\U95;
figure
plot(0:1e3, sqrt(conj(U95E).*U95E))
plot(0:1e3, log10(sqrt(conj(U95E).*U95E)))
[U94, mniter94, matvecs94, est_errors94, history94] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 5e3, 9, 4, 1e-5, 10, 16, false);
est_errors94
norm(U94(:, end)-Uex(:,2))/norm(Uex(:,2))
5.6923e-04/1.2870e-05
est_errors95
norm(U95(:, end)-Uex(:,2))/norm(Uex(:,2))
9.9553e-05/1.9746e-06
[U93, mniter93, matvecs93, est_errors93, history93] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 5e3, 9, 3, 1e-5, 10, 16, false);
est_errors93
norm(U93(:, end)-Uex(:,2))/norm(Uex(:,2))
3.8802e-03/4.0928e-05
[U3e3, mniter3e3, matvecs3e3, est_errors3e3, history3e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 3e3, 9, 7, 1e-5, 10, 16, false);
est_errors3e3
norm(U3e3(:, end)-Uex(:,2))/norm(Uex(:,2))
2.5494e-03/8.4387e-05
mniter3e3
matvecs3e3
matvecs
[U3e3_95, mniter3e3_95, matvecs3e3_95, est_errors3e3_95, history3e3_95] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 3e3, 9, 5, 1e-5, 10, 16, false);
norm(U3e3_95(:, end)-Uex(:,2))/norm(Uex(:,2))
est_errors3e3_95
figure
plot(0:1e3, vecnorm(U3e3_95))
viewP(U3e3_59, x240, 0.01)
viewP(U3e3_95, x240, 0.01)
viewP(U3e3, x240, 0.01)
[U3e3_96, mniter3e3_96, matvecs3e3_96, est_errors3e3_96, history3e3_96] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 3e3, 9, 6, 1e-5, 10, 16, false);
est_errors3e3_96
norm(U3e3_96(:, end)-Uex(:,2))/norm(Uex(:,2))
2.3509e-02/4.4931e-04
viewP(U3e3_96, x240, 0.01)
[U1e3, mniter1e3, matvecs1e3, est_errors1e3, history1e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 7, 1e-5, 10, 16, false);
[U1e3, mniter1e3, matvecs1e3, est_errors1e3, history1e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 7, 1e-5, 10, 16);
[U1e3, mniter1e3, matvecs1e3, est_errors1e3, history1e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 9, 1e-5, 10, 16, false);
[U1e3, mniter1e3, matvecs1e3, est_errors1e3, history1e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 11, 1e-5, 10, 16, false);
clear U1e3 mniter1e3 matvecs1e3 est_errors1e3 history1e3
[U1e3, mniter1e3, matvecs1e3, est_errors1e3, history1e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 13, 1e-5, 10, 16, false);
[U2e3, mniter2e3, matvecs2e3, est_errors2e3, history2e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 7, 1e-5, 10, 16, false);
[U2e3, mniter2e3, matvecs2e3, est_errors2e3, history2e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 1e3, 9, 9, 1e-5, 10, 16, false);
[U2e3, mniter2e3, matvecs2e3, est_errors2e3, history2e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 2e3, 9, 7, 1e-5, 10, 16, false);
[U2e3, mniter2e3, matvecs2e3, est_errors2e3, history2e3] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 2e3, 9, 9, 1e-5, 10, 16, false);
norm(U2e3(:, end)-Uex(:,2))/norm(Uex(:,2))
viewP(U2e3, x240, 0.01)
[U2e3_11, mniter2e3_11, matvecs2e3_11, est_errors2e3_11, history2e3_11] = SemiGlobal1(Gop, Gdiff_op, 0, [], [], fi0240, 0:1e3, 2e3, 9, 11, 1e-5, 10, 16, false);
norm(U2e3_11(:, end)-Uex(:,2))/norm(Uex(:,2))
clear all