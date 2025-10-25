test_harmonic
[U, mniter, matvecs, max_errors, history] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - ones(1, Nt_ts)*cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
[Ud, mniterd, matvecsd, max_errorsd, historyd] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - ones(1, Nt_ts)*cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
t(s_ext_i)
[U1, mniter1, matvecs1, max_errors1, history1] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
max(max(abs(U-U1)))
clear U1 mniter1 matvecs1 max_errors1 history1
[Ud, mniterd, matvecsd, max_errorsd, historyd] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
max(max(abs(U-Ud)))
vecnorm(Ud)
historyd.niter
[Ud, mniterd, matvecsd, max_errorsd, historyd] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
t_2ts(1:total_Nt_ts)
[Ud, mniterd, matvecsd, max_errorsd, historyd] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
max(max(abs(U-Ud)))
max_errors
max_errorsd
[Ud, mniterd, matvecsd, max_errorsd, historyd] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
max_errorsd
historyd.texp_error
historyd.texp_error.safe
historyd.texp_error.safe1
historyd.texp_error.safe1 - historyd.texp_error.safe
abs(historyd.texp_error.safe1 - historyd.texp_error.safe)./abs(historyd.texp_error.safe)
abs(historyd.texp_error.safe)
[Ud, mniterd, matvecsd, max_errorsd, historyd] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
historyd.texp_error.exact
historyd.texp_error.exact./historyd.texp_error.safe
[U, mniter, matvecs, errors, Ferror, texp_dif] = TDHxp_tests(K, V, @(u, x, t) x*cos(t), [], [0 195], fi0, x, [0 T], Nts, Nt_ts, Ncheb, tol);
[U, mniter, matvecs, max_errors, history] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
[Uo, mnitero, matvecso, errorso, Ferroro, texp_difo] = TDHxp_tests(K, V, @(u, x, t) x*cos(t), [], [0 195], fi0, x, [0 T], Nts, Nt_ts, Ncheb, tol);
max(max(abs(U-Uo)))
size(Uo)
max(max(abs(U-Uo(:, 1:2:end))))
errorso(11,:)
errorso(:, 11)
historyd.texp_error.exact./errors0(:,11)
historyd.texp_error.exact./errorso(:,11)
historyd.texp_error.exact./errorso(:,11).'
[historyd.texp_error.exact.'./errorso(:,11) historyd.texp_error.exact.' errorso(:,11)]
(historyd.texp_error.exact + 1e-18)./errorso(:,11).'
[U5d, mniter5d, matvecs5d, max_errors5d, history5d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 5, 13, 1e-10);
[U5o, mniter5o, matvecs5o, errors5o, Ferror5o, texp_dif5o] = TDHxp_tests(K, V, @(u, x, t) x*cos(t), [], [0 195], fi0, x, [0 T], Nts, 5, 13, 1e-10);
max_errors5d
history5d.texp_error.exact
sum(history5d.texp_error.exact)
errors5o(:, 11).'
history5d.texp_error.exact./errors5o(:,11).'
mx5 = conj(U).*U.*x;
mx5 = sum(conj(U).*U.*x);
max(abs(mx5-mx))
mx5 = sum(conj(U5d).*U5d.*x);
max(abs(mx5-mx))
sum(history5d.texp_error.exact)
max_errors5d
historyd5.texp_error.safe1-historyd5.texp_error.safe
history5d.texp_error.safe1-historyd5.texp_error.safe
history5d.texp_error.safe1-history5d.texp_error.safe
(history5d.texp_error.safe1-history5d.texp_error.safe)./history5d.texp_error.safe1
[U6o, mniter6o, matvecs6o, errors6o, Ferror6o, texp_dif6o] = TDHxp_tests(K, V, @(u, x, t) x*cos(t), [], [0 195], fi0, x, [0 T], Nts, 6, 13, 1e-10);
[U6d, mniter6d, matvecs6d, max_errors6d, history6d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 6, 13, 1e-10);
max_errors6d
max_errors5d
sum(history6d.texp_error.exact)
mx6 = sum(conj(U6d).*U6d.*x);
max(abs(mx5-error))
max(abs(mx5 - (-0.5*sin(t).*t)))
max(abs(mx6 - (-0.5*sin(t).*t)))
history6d.texp_error.exact - errors6o(:, 7).'
(history6d.texp_error.exact - errors6o(:, 7).')./history6d.texp_error.exact
mx_ex = (-0.5*sin(t).*t);
mp_ex = -0.5*(sin(t) + t.*cos(t));
mx_new = evmu(U, x);
max(abs(mx-mx_new))
max(imag(mx-mx_new))
max(imag(mx_new))
max(imag(mx))
mp = evp(U, p);
figure
whos
plot(t, mx)
mx = evmu(U, x);
plot(t, mx)
hold on
plot(t, mp)
max(abs(mp-mp_ex))
max(abs(mx-mx_ex))
psi_exact = pi^(-1/4)*exp(1i*(mp_exact*(x - mx_exact/2) - t/2) - (x - mx_exact).^2/2)*sqrt(dx);
psi_exact = pi^(-1/4)*exp(1i*(mp_ex*(x - mx_ex/2) - t/2) - (x - mx_ex).^2/2)*sqrt(dx);
psi_exact = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex/2) - t/2) - (x - mx_ex).^2/2)*sqrt(dx);
size(psi_exact)
clear psi_exact
Uex = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex/2) - t/2) - (x - mx_ex).^2/2)*sqrt(dx);
max(vecnorm(Uex-U))
vecnorm(Uex).^2
Uex = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex) - t/2) - (x - mx_ex).^2/2)*sqrt(dx);
max(vecnorm(Uex-U))
Uex = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex/2) - t/2) - (x - mx_ex).^2/2)*sqrt(dx);
max(abs(mx_ex - evmu(Uex, x)))
max(abs(mp_ex - evp(Uex, p)))
Uex = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex/2)) - (x - mx_ex).^2/2)*sqrt(dx);
max(vecnorm(Uex-U))
Uex = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex)) - (x - mx_ex).^2/2)*sqrt(dx);
max(vecnorm(Uex-U))
Uex = pi^(-1/4)*exp(1i*(mp_ex.*(x - mx_ex/2)) - (x - mx_ex).^2/2)*sqrt(dx);
phaseU = U./Uex
angle(phaseUex(65,:))
angle(phaseU(65,:))
angleU = -angle(phaseU(65,:));
figure
plot(t, angleU)
0.7*83+30
x = 53/0.7
angleUfixed = phi_cont1(angleU);
figure
plot(t, angleUfixed)
Uex1 = pi^(-1/4)*exp(1i*(mp_ex.*x) - (x - mx_ex).^2/2)*sqrt(dx);
phaseU1 = U./Uex1
vecnorm(Uex1)
Uex1 = pi^(-1/4)*exp(1i*(x.*mp_ex) - (x - mx_ex).^2/2)*sqrt(dx);
vecnorm(Uex1)
x.*mp_ex
Uex1 = pi^(-1/4)*exp(1i*(x*mp_ex) - (x - mx_ex).^2/2)*sqrt(dx);
vecnorm(Uex1)
x*mp_ex
size(x)
test_harmonic
Uex1 = pi^(-1/4)*exp(1i*(x*mp_ex) - (x - mx_ex).^2/2)*sqrt(dx);
vecnorm(Uex1)
phaseU1 = U./Uex1;
angleU1 = -angle(phaseU(65,:));
figure
plot(t, angleU1)
angleU1 = -angle(phaseU1(65,:));
plot(t, angleU1)
angleU1fixed = phi_cont1(angleU1);
plot(t, angleU1fixed)
[Uex_num, mniterex, matvecsex, max_errorsex, historyex] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 2e3, 9, 13, eps);
1-vecnorm(Uex_num)
mniterex
matvecsex
max_errorsex
1-vecnorm(U)
sum(history5d.texp_error.exact)
sum(history6d.texp_error.exact)
[U5d, mniter5d, matvecs5d, max_errors5d, history5d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 2e3, 5, 13, eps);
[U5d, mniter5d, matvecs5d, max_errors5d, history5d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 2e3, 5, 13, eps, 10, 16, false);
[U5d, mniter5d, matvecs5d, max_errors5d, history5d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 5, 13, eps, 10, 16, false);
sum(history5d.texp_error.exact)
[U6d, mniter6d, matvecs6d, max_errors6d, history6d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 6, 13, eps);
[U6d, mniter6d, matvecs6d, max_errors6d, history6d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 6, 13, eps, 10, 16, false);
mniter6d
mniter5d
sum(history6d.texp_error.exact)
vecnorm(U5d-Uex_num)
sum(history6d.texp_error.exact)
sum(history5d.texp_error.exact)
vecnorm(U5d-Uex_num) - history5d.texp_error.exact
texp_er5num = vecnorm(U5d-Uex_num);
texp_er5num - history5d.texp_error.exact
texp_er5num(2:end) - history5d.texp_error.exact
size(texp_er5num)
size(history5d.texp_error.exact)
texp_er5num(2:end) - history5d.texp_error.exact(2:2:end)
cum_texper5 = cumsum(history5d.texp_error.exact);
texp_er5num(2:end) - cum_texper5(2:2:end)
size(cum_texper5)
(cum_texper5(2:2:end) - texp_er5num(2:end))./texp_er5num(2:end)
figure
plot(t(2:end), (cum_texper5(2:2:end) - texp_er5num(2:end))./texp_er5num(2:end))
figure
plot(t, texp_er5num)
plot(t, log10(texp_er5num))
hold on
mx5 = sum(conj(U5d).*U5d.*x);
mx6 = sum(conj(U6d).*U6d.*x);
plot(t, log10(abs(mx5-mx_ex)))
mp5 = evp(U5d, p);
plot(t, log10(abs(mp5-mp_ex)))
texp_er5x = vecnorm(abs(U5d)-abs(Uex));
plot(t, log10(texp_er5x))
figure
plot(t, mx_ex)
save errors_harmonic
texp_er6num = vecnorm(U6d-Uex_num);
cum_texper6 = cumsum(history6d.texp_error.exact);
figure
plot(t(2:end), (cum_texper6(2:2:end) - texp_er6num(2:end))./texp_er6num(2:end))
[U7d, mniter7d, matvecs7d, max_errors7d, history7d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 7, 13, eps, 10, 16, false);
sum(history7d.texp_error.exact)
texp_er7num = vecnorm(U7d-Uex_num);
cum_texper7 = cumsum(history7d.texp_error.exact);
figure
plot(t(2:end), (cum_texper7(2:2:end) - texp_er7num(2:end))./texp_er7num(2:end))
cum_texper7
[cum_texper7(2:2:end).', texp_er7num(2:end).']
[U8d, mniter8d, matvecs8d, max_errors8d, history8d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 8, 13, eps, 10, 16, false);
sum(history8d.texp_error.exact)
texp_er8num = vecnorm(U8d-Uex_num);
cum_texper8 = cumsum(history8d.texp_error.exact);
figure
plot(t(2:end), (cum_texper8(2:2:end) - texp_er8num(2:end))./texp_er8num(2:end))
[cum_texper8(2:2:end).', texp_er8num(2:end).']
nnz(cum_texper8>1e-15)
nnz(cum_texper8>1e-14)
save errors_harmonic
whos
angle_analytic = -(t/2 - (sin(2*t)/4 - t.*cos(2*t))/8);
figure
plot(t, angle_analytic)
angle_analytic = (t/2 - (sin(2*t)/4 - t.*cos(2*t))/8);
plot(t, angle_analytic)
angle_analytic = (t/2 - (sin(2*t)/4 - t.*cos(2*t)/2)/8);
plot(t, angle_analytic)
hold on
plot(t, angle_analytic)
clear historyex.U
figure
plot(t, abs(angleUfixed - angle_analytic))
plot(t, log10(abs(angleUfixed - angle_analytic)))
figure
plot(t, log10(abs(exp(1i*abs(angleUfixed - angle_analytic))-1)))
save errors_harmonic
Uex_phase = pi^(-1/4)*exp(-1i*angle_analytic).*exp(1i*(mp_ex.*(x - mx_ex/2)) - (x - mx_ex).^2/2)*sqrt(dx);
max(abs(vecnorm(Uex_num - Uex_phase)))
texp_er8anal = vecnorm(U8d-Uex_phase);
figure
plot(t(2:end), (cum_texper8(2:2:end) - texp_er8anal(2:end))./texp_er8anal(2:end))
hold on
plot(t(2:end), (cum_texper8(2:2:end) - texp_er8num(2:end))./texp_er8num(2:end))
[cum_texper8(2:2:end).', texp_er8num(2:end).' texp_er8anal(2:end).']
[t(2:end).' cum_texper8(2:2:end).', texp_er8num(2:end).' texp_er8anal(2:end).']
max((vecnorm(Uex_num - Uex_phase))
max((vecnorm(Uex_num - Uex_phase)))
figure
plot(t, vecnorm(Uex_num - Uex_phase))
figure
plot(t, vecnorm(Uex_num - U8d))
plot(t(1:101), vecnorm(Uex_num(:,1:101) - U8d(1:101)))
plot(t(1:101), vecnorm(Uex_num(:,1:101) - U8d(:,1:101)))
plot(t(1:51), vecnorm(Uex_num(:,1:51) - U8d(:,1:51)))
plot(t(1:51), vecnorm(Uex_phase(:,1:51) - U8d(:,1:51)))
hold on
plot(t(1:51), vecnorm(Uex_phase(:,1:51) - Uex_num(:,1:51)))
hold on
plot(t(1:51), vecnorm(Uex_phase(:,1:51) - U8d(:,1:51)))
figure
plot(t, log10(vecnorm(Uex_phase - U5d)))
hold on
plot(t, log10(vecnorm(Uex_num - U5d)))
[U5d2, mniter5d2, matvecs5d2, max_errors5d2, historyd2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 5, 17, eps, 10, 16, false);
[U5d2, mniter5d2, matvecs5d2, max_errors5d2, historyd2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 5, 30, eps, 10, 16, false);
mniter5d2
max_errors5d2
figure
plot(t, log10(vecnorm(Uex_num - U5d2)))
sum(history5d2.texp_error.exact)
[U5d, mniter5d, matvecs5d, max_errors5d, history5d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 5, 13, eps, 10, 16, false);
[U5d2, mniter5d2, matvecs5d2, max_errors5d2, history5d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 5, 30, eps, 10, 16, false);
clear historyd2
sum(history5d2.texp_error.exact)
cum_texper5_2 = cumsum(history5d2.texp_error.exact);
hold on
texp_er5_2anal = vecnorm(U5d2-Uex_phase);
clf
plot(t, tex_er5_2anal)
plot(t, texp_er5_2anal)
plot(t, log10(texp_er5_2anal))
hold on
plot(t, log10(cum_texper5_2))
plot(t(2:end), log10(cum_texper5_2))
figure
plot(t(2:end), cum_texper5_2-texp_er5_2anal(2:end))
plot(t(2:end), (cum_texper5_2-texp_er5_2anal(2:end))./texp_er5_2anal(2:end))
[U5d3, mniter5d3, matvecs5d3, max_errors5d3, history5d3] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 50, 5, 30, eps, 10, 16, false);
cum_texper5_2 = cumsum(history5d3.texp_error.exact);
texp_er5_3anal = vecnorm(U5d3-Uex_phase);
figure
plot(t(2:end), (cum_texper5_3-texp_er5_3anal(2:end))./texp_er5_3anal(2:end))
cum_texper5_2 = cumsum(history5d2.texp_error.exact);
cum_texper5_3 = cumsum(history5d3.texp_error.exact);
plot(t(2:end), (cum_texper5_3-texp_er5_3anal(2:end))./texp_er5_3anal(2:end))
plot(t(4:2:end), (cum_texper5_3-texp_er5_3anal(2:end))./texp_er5_3anal(2:end))
plot(0.2:0.2:10, (cum_texper5_3-texp_er5_3anal(2:end))./texp_er5_3anal(2:end))
t50 = 0:0.2:10;
size(cum_texper5)
size(cum_texper5_3)
plot(t50(2:end), (cum_texper5_3-texp_er5_3anal(3:2:end))./texp_er5_3anal(3:2:end))
sum(history5d3.texp_error.exact)
[U6d2, mniter6d2, matvecs6d2, max_errors6d2, history6d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 5, 30, eps, 10, 16, false);
cum_texper6_2 = cumsum(history6d2.texp_error.exact);
cum_texper6_2(end)
[U6d2, mniter6d2, matvecs6d2, max_errors6d2, history6d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 6, 30, eps, 10, 16, false);
cum_texper6_2 = cumsum(history6d2.texp_error.exact);
cum_texper6_2(end)
figure
texp_er6_2anal = vecnorm(U6d2-Uex_phase);
plot(t(2:end), (cum_texper6_2-texp_er6_2anal(2:end))./texp_er6_2anal(2:end))
figure
plot(t(2:end), cum_texper6)
plot(t(2:end), cum_texper6_2)
hold on
plot(t, texp_er6_2anal)
[U7d2, mniter7d2, matvecs7d2, max_errors7d2, history7d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 7, 30, eps, 10, 16, false);
cum_texper7_2 = cumsum(history7d2.texp_error.exact);
cum_texper7_2(end)
texp_er7_2anal = vecnorm(U7d2-Uex_phase);
figure
plot(t(2:end), (cum_texper7_2-texp_er7_2anal(2:end))./texp_er7_2anal(2:end))
[U8d2, mniter8d2, matvecs8d2, max_errors8d2, history8d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 8, 30, eps, 10, 16, false);
cum_texper8_2 = cumsum(history8d2.texp_error.exact);
figure
plot(t(2:end), (cum_texper8_2-texp_er8_2anal(2:end))./texp_er8_2anal(2:end))
texp_er8_2anal = vecnorm(U8d2-Uex_phase);
plot(t(2:end), (cum_texper8_2-texp_er8_2anal(2:end))./texp_er8_2anal(2:end))
[U6d3, mniter6d3, matvecs6d3, max_errors6d3, history6d3] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 50, 6, 30, eps, 10, 16, false);
mniter6d3
mniter5d3
sum(history6d3.texp_error.exact)
figure
texp_er6_3anal = vecnorm(U6d3-Uex_phase);
cum_texper6_3 = cumsum(history6d3.texp_error.exact);
plot(t(2:end), (cum_texper6_3-texp_er6_3anal(2:end))./texp_er6_3anal(2:end))
plot(t50(2:end), (cum_texper6_3-texp_er6_3anal(3:2:end))./texp_er6_3anal(3:2:end))
[~, ~, ~, E, P, H] = gsV(@(x) x.^2/2, [-L/2, L/2], 128);
E
E = real(E);
Uex_phaseE = P\Uex_phase;
figure
plot(0:127, Uex_phaseE(:,end).*conj(Uex_phaseE(:,end)))
E(30)
30*0.2
30*0.05
sum(Uex_phaseE(:,end).*conj(Uex_phaseE(:,end)).*E)
20*0.05
14*0.05
200*0.05
%-- 08/07/2022 15:04 --%
load errors_harmonic
Uex_phase = pi^(-1/4)*exp(-1i*angle_analytic).*exp(1i*(mp_ex.*(x - mx_ex/2)) - (x - mx_ex).^2/2)*sqrt(dx);
max(abs(vecnorm(Uex_num - Uex_phase)))
texp_er8anal = vecnorm(U8d-Uex_phase);
[U5d2, mniter5d2, matvecs5d2, max_errors5d2, history5d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 5, 30, eps, 10, 16, false);
cum_texper5_2 = cumsum(history5d2.texp_error.exact);
texp_er5_2anal = vecnorm(U5d2-Uex_phase);
texp_er5_3anal = vecnorm(U5d3-Uex_phase);
[U5d3, mniter5d3, matvecs5d3, max_errors5d3, history5d3] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 50, 5, 30, eps, 10, 16, false);
cum_texper5_3 = cumsum(history5d3.texp_error.exact);
texp_er5_3anal = vecnorm(U5d3-Uex_phase);
t50 = 0:0.2:10;
[U6d2, mniter6d2, matvecs6d2, max_errors6d2, history6d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 6, 30, eps, 10, 16, false);
cum_texper6_2 = cumsum(history6d2.texp_error.exact);
texp_er6_2anal = vecnorm(U6d2-Uex_phase);
[U7d2, mniter7d2, matvecs7d2, max_errors7d2, history7d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 7, 30, eps, 10, 16, false);
cum_texper7_2 = cumsum(history7d2.texp_error.exact);
cum_texper7_2(end)
texp_er7_2anal = vecnorm(U7d2-Uex_phase);
[U8d2, mniter8d2, matvecs8d2, max_errors8d2, history8d2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 100, 8, 30, eps, 10, 16, false);
cum_texper8_2 = cumsum(history8d2.texp_error.exact);
texp_er8_2anal = vecnorm(U8d2-Uex_phase);
[U6d3, mniter6d3, matvecs6d3, max_errors6d3, history6d3] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, 50, 6, 30, eps, 10, 16, false);
texp_er6_3anal = vecnorm(U6d3-Uex_phase);
cum_texper6_3 = cumsum(history6d3.texp_error.exact);
[~, ~, ~, E, P, H] = gsV(@(x) x.^2/2, [-L/2, L/2], 128);
E = real(E);
Uex_phaseE = P\Uex_phase;
save errors_harmonic
figure
%-- 11/07/2022 14:59 --%
load errors_harmonic
[UK5d, mniterK5d, matvecsK5d, max_errorsK5d, historyK5d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 5, eps, 10, 16, false);
norm(UK5d(:, end))
1-norm(UK5d(:, end))
sum(historyK5d.texp_error.exact)
sum(historyK5d.fUerror)
fUer5anal = vecnorm(UK5d-Uex_phase);
figure
plot(t, fUer5anal)
historyK5d.fUerror
40^9/factorial(9)
historyK5d.f_error
UK5dE = P\UK5d;
figure
plot(t, UK5dE.*conj(UK5dE))
plot(t, log10(UK5dE.*conj(UK5dE)))
figure
cum_fUerK5d = cumsum(historyK5d.fUerror);
hold on
plot(t, cum_fUerK5d)
plot(t, cum_fUerK5d(1:2:end))
plot(t(2:end), cum_fUerK5d(1:2:end))
hold on
plot(t(2:end), cum_fUerK5d(1:2:end))
plot(t(2:end), cum_fUerK5d(1:2:end)-fUerK5d_anal(2:end))
plot(t(2:end), cum_fUerK5d(1:2:end)-fUer5anal(2:end))
plot(t(2:end), (cum_fUerK5d(1:2:end)-fUer5anal(2:end))./fUer5anal)
plot(t(2:end), (cum_fUerK5d(1:2:end)-fUer5anal(2:end))./fUer5anal(2:end))
[t(2:end).', cum_fUerK5d(1:2:end), fUer5anal(2:end)]
[t(2:end).', cum_fUerK5d(1:2:end).', fUer5anal(2:end).']
[UK5d_2, mniterK5d_2, matvecsK5d_2, max_errorsK5d_2, historyK5d_2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 7, 5, eps, 10, 16, false);
norm(UK5d(:, end))
sum(historyK5d_2.fUerror)
norm(UK5d_2(:, end))
[UK5d_2, mniterK5d_2, matvecsK5d_2, max_errorsK5d_2, historyK5d_2] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 7, 7, eps, 10, 16, false);
norm(UK5d_2(:, end))
sum(historyK5d.fUerror)
sum(historyK5d_2.fUerror)
sum(historyK5d_2.texp_error.exact)
40^7/factorial(7)
clear UK5d_2 mniterK5d_2 matvecsK5d_2 max_errorsK5d_2 historyK5d_2
[
]
[UK7d, mniterK7d, matvecsK7d, max_errorsK7d, historyK7d] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 7, 7, eps, 10, 16, false);
cum_fUerK7d = cumsum(historyK7d.fUerror);
fUer7anal = vecnorm(UK7d-Uex_phase);
UK7dE = P\UK7d;
figure
plot(t, log10(UK7dE.*conj(UK7dE)))
figure
plot(t(2:end), (cum_fUerK7d(1:2:end)-fUer7anal(2:end))./fUer7anal(2:end))
save errors_harmonic
[UK5, mniterK5, matvecsK5, max_errorsK5, historyK5] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 5, eps, 10, 16, false);
sum(historyK5d.fUerror)
sum(historyK5.fUerror)
figure
cum_fUerK5 = cumsum(historyK5.fUerror);
plot(t(2:end), (cum_fUerK5d(1:2:end)-fUer5anal(2:end))./fUer5anal(2:end))
hold on
plot(t(2:end), (cum_fUerK5(1:2:end)-fUer5anal(2:end))./fUer5anal(2:end))
historyK5.f_error
max_errorsK5
historyK5d.f_error
figure
plot(t, historyK5d.f_error)
plot(t(2:end), historyK5d.f_error)
plot(0.05:0.05:10, historyK5d.f_error)
[UK5ar, mniterK5ar, matvecsK5ar, max_errorsK5ar, historyK5ar] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, Nts, 9, 5, eps, 10, 16, false);
UK5dEar = P\UK5ar;
figure
plot(t, log10(UK5arE.*conj(UK5arE)))
clear UK5dEar
UK5arE = P\UK5ar;
plot(t, log10(UK5arE.*conj(UK5arE)))
figure
cum_fUerK5ar = cumsum(historyK5ar.fUerror);
fUer5ar_anal = vecnorm(UK5ar-Uex_phase);
plot(t, fUer5anal)
hold on
plot(t, fUer5ar_anal)
figure
plot(t(2:end), (cum_fUerK5ar(1:2:end)-fUer5anal(2:end))./fUer5anal(2:end))
plot(t(2:end), (cum_fUerK5ar(1:2:end)-fUer5ar_anal(2:end))./fUer5ar_anal(2:end))
[t(2:end).', cum_fUerK5ar(1:2:end).', fUer5ar_anal(2:end).']
[UK7ar, mniterK7ar, matvecsK7ar, max_errorsK7ar, historyK7ar] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, t, Nts, 7, 7, eps, 10, 16, false);
fUer7ar_anal = vecnorm(UK7ar-Uex_phase);
cum_fUerK7ar = cumsum(historyK7ar.fUerror);
figure
plot(t(2:end), (cum_fUerK7ar(1:2:end)-fUer7ar_anal(2:end))./fUer7ar_anal(2:end))
[t(2:end).', cum_fUerK7ar(1:2:end).', fUer7ar_anal(2:end).']
sum(historyK7ar.texp_error.exact)
clear UK7ar mniterK7ar matvecsK7ar max_errorsK7ar historyK7ar
clear  fUer7ar_anal cum_fUerK7ar
save errors_harmonic
[Ud1i, mniterd1i, matvecs1i, max_errors1i, history1i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [-195*1i, 0], [], fi0, t, Nts, 9, 13, eps, 10, 1, false);
[Ud1i, mniterd1i, matvecs1i, max_errors1i, history1i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, eps, 10, 1, false);
norm(Ud1i(:, end))
conv_er1i = vecnorm(Ud1i-Uex_phase);
figure
plot(t, conv_er1i)
figure
plot(t, conv_er1i)
plot(t(2:end), conv_er1i(2:end))
[Ud1i, mniterd1i, matvecs1i, max_errors1i, history1i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, eps, 1, 16, false);
norm(Ud1i(:, end))
conv_er1i = vecnorm(Ud1i-Uex_phase);
figure
plot(t, conv_er1i)
sum(historyd1i.conv_error.safe)
clear matvecs1i max_errors1i history1i
[Ud1i, mniterd1i, matvecsd1i, max_errorsd1i, historyd1i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, eps, 1, 16, false);
sum(historyd1i.conv_error.safe)
sum(historyd1i.conv_error.exact)
cum_conv_er_ex_d1i = cumsum(historyd1i.conv_error.exact);
hold on
plot(t, cum_conv_er_ex_d1i)
plot(0.05:0.05:10, cum_conv_er_ex_d1i)
figure
plot(t(2:end), (cum_conv_er_ex_d1i(1:2:end)-conv_er1i(2:end))./conv_er1i(2:end))
[Udc, mniterdc, matvecsdc, max_errorsdc, historydc] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, 1e-5, 10, 16, false);
[Udc, mniterdc, matvecsdc, max_errorsdc, historydc] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, 1e-8, 10, 16, false);
sum(historydc.conv_error.exact)
sum(historydc.texp_error.exact)
sum(historydc.fUerror)
figure
conv_erdc = vecnorm(Udc-Uex_phase);
plot(t, conv_erdc)
hold on
cum_conv_er_ex_dc = cumsum(historydc.conv_error.exact);
hold on
plot(0.05:0.05:10, conv_erdc)
plot(0.05:0.05:10, cum_conv_er_ex_dc)
[Uco, mniterco, matvecsco, errorsco, Ferrorco, texp_difco] = TDHxp_tests(K, V, @(u, x, t) x*cos(t), [], [0 195], fi0, x, [0 T], Nts, 9, 13, 1e-8);
matvecsdc
matvecsco
mniterco
mniterdc
[U5d_temp, mniter5d_temp, matvecs5d_temp, max_errors5d_temp, history5d_temp] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 5, 13, eps, 10, 16, false);
mniter5d_temp
mniter5d
[U5d_temp, mniter5d_temp, matvecs5d_temp, max_errors5d_temp, history5d_temp] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 5, 13, eps, 10, 16, false);
mniter5d_temp
mniter5d
[Udc, mniterdc, matvecsdc, max_errorsdc, historydc] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, 1e-8, 10, 16, false);
mniterdc
cum_conv_er_ex_dc = cumsum(historydc.conv_error.exact);
save errors_harmonic
max(max(abs(U5d-U5d_temp)))
size(max_errorsdc)
size(errorsco)
max(max(abs(Udc-Uco)))
max(max(abs(Udc-Uco(:, 1:2:end))))
errorsco(4,:) - historydc.conv_error.exact
errorsco(4,:).' - historydc.conv_error.exact
size(historydc.conv_error.exact)
errorsco(:,4).' - historydc.conv_error.exact
errorsco(:,3).' - historydc.conv_error.safe
conv_erdc = vecnorm(Udc-Uex_phase);
clf
plot(t, conv_erdc)
hold on
cum_conv_er_ex_dc = cumsum(historydc.conv_error.exact);
plot(0.05:0.05:10, cum_conv_er_ex_dc)
figure
plot(t, conv_erdc)
hold on
cum_conv_er_s_dc = cumsum(historydc.conv_error.safe);
plot(0.05:0.05:10, cum_conv_er_s_dc)
figure
plot(0.05:0.05:10, cum_conv_er_s_dc./cum_conv_er_ex_dc)
figure
plot(t(2:end), (cum_texper5(2:2:end) - texp_er5num(2:end))./texp_er5num(2:end))
save errors_harmonic
[Ud1i, mniterd1i, matvecsd1i, max_errorsd1i, historyd1i] = SemiGlobalDirty(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, 9, 13, eps, 1, 16, false);
cum_conv_er_ex_d1i = cumsum(historyd1i.conv_error.exact);
cum_conv_er_s_d1i = cumsum(historyd1i.conv_error.safe);
conv_er1i = vecnorm(Ud1i-Uex_phase);
figure
plot(t, conv_er1i)
hold on
plot(0.05:0.05:10, cum_conv_er_s_d1i)
plot(0.05:0.05:10, cum_conv_er_ex_d1i)
figure
plot(t(2:end), (cum_conv_er_ex_d1i(1:2:end)-conv_er1i(2:end))./conv_er1i(2:end))
[t(2:end).', cum_conv_er_ex_d1i(1:2:end).', conv_er1i(2:end).']
figure
plot(t(2:end), (cum_conv_er_s_d1i(1:2:end)-conv_er1i(2:end))./conv_er1i(2:end))
figure
plot(t(2:end), cum_conv_er_s_d1i(1:2:end)./cum_conv_er_ex_d1i(1:2:end))
save errors_harmonic
%%%%%
load errors_harmonic
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase, [0 10], 9, 9, 150, 13, 1);
150*10^0.1
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase, [0 10], 9, 9, 200, 13, 1);
size(Uex_phase)
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 9, 200, 13, 1);
hold on
plot(log10(allmv99), log10(est_ers99.texp_exact), '-o')
plot(log10(allmv99), log10(est_ers99.fU), '-o')
plot(log10(allmv99), log10(est_ers99.conv), '-o')
plot(log10(allmv99), log10(est_ers99.conv_exact), '-o')
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 7, 9, 200, 10, 1);
hold on
plot(log10(allmv79), log10(est_ers79.texp_exact), '-o')
plot(log10(allmv99), log10(est_ers99.texp_exact), '-o')
plot(log10(allmv99), log10(est_ers99.fU), '-o')
plot(log10(allmv99), log10(est_ers99.conv_exact), '-o')
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 9, 200, 7, 1);
[allNt79, allmv79, aller79, est_ers79] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 7, 9, 200, 10, 1);
[allNt7i2, allmv7i2, aller7i2, est_ers7i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 7, 13, 200, 13, 2);
hold on
plot(log10(allmv7i2), log10(est_ers7i2.texp_exact), '-o')
plot(log10(allmv7i2), log10(est_ers7i2.fU), '-o')
plot(log10(allmv7i2), log10(est_ers7i2.conv_exact), '-o')
[allNt5i2, allmv5i2, aller5i2, est_ers5i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 5, 13, 200, 10, 2);
[allNt5i2, allmv5i2, aller5i2, est_ers5i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 5, 13, 200, 13, 2);
hold on
plot(log10(allmv5i2), log10(est_ers5i2.texp_exact), '-o')
plot(log10(allmv5i2), log10(est_ers5i2.fU), '-o')
plot(log10(allmv5i2), log10(est_ers5i2.conv_exact), '-o')
[allNt6i2, allmv6i2, aller6i2, est_ers6i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 6, 13, 200, 13, 2);
hold on
plot(log10(allmv6i2), log10(est_ers6i2.texp_exact), '-o')
plot(log10(allmv6i2), log10(est_ers6i2.fU), '-o')
plot(log10(allmv6i2), log10(est_ers6i2.conv_exact), '-o')
polyfit(log10(allmv5i2(1:12)), log10(aller5i2.texp_exact(1:12)), 1)
polyfit(log10(allmv5i2(1:12)), log10(aller5i2(1:12)), 1)
polyfit(log10(allmv6i2(1:9)), log10(aller6i2(1:9)), 1)
polyfit(log10(allmv6i2(2:9)), log10(aller6i2(2:9)), 1)
save errors_harmonic
figure
plot(log10(allmv79), log10(aller79), '-o')
[allNt55, allmv55, aller55, est_ers55] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 5, 5, 200, 16, 1);
hold on
plot(log10(allmv55), log10(est_ers55.texp_exact), '-o')
plot(log10(allmv55), log10(est_ers55.fU), '-o')
plot(log10(allmv55), log10(est_ers55.conv_exact), '-o')
polyfit(log10(allmv6i2(2:9)), log10(aller6i2(2:9)), 1)
polyfit(log10(allmv55(1:12)), log10(aller55(1:12)), 1)
polyfit(log10(allmv55(2:12)), log10(aller55(2:12)), 1)
polyfit(log10(allmv55(2:11)), log10(aller55(2:11)), 1)
polyfit(log10(allmv5i2(1:12)), log10(aller5i2(1:12)), 1)
polyfit(log10(allmv79(1:12)), log10(aller79(1:12)), 1)
polyfit(log10(allmv79(1:7)), log10(aller79(1:7)), 1)
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 5, 200, 7, 1);
[allNt99, allmv99, aller99, est_ers99] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 9, 200, 7, 1);
hold on
plot(log10(allmv99), log10(est_ers99.conv_exact), '-o')
polyfit(log10(allmv99(1:5)), log10(aller99(1:5)), 1)
plot(log10(allmv99), log10(est_ers99.fU), '-o')
plot(log10(allmv99), log10(est_ers99.conv_exact), '-o')
plot(log10(allmv99), log10(est_ers99.texp_exact), '-o')
save errors_harmonic
%%%%%%%%%%%%%%%
load conv_ratios
load errors_harmonic
whos
[est_ers79.conv_cheap*2*all_gamma(1,6), est_ers79.conv_exact, aller79]
[est_ers79.conv_cheap(1:8)*2*all_gamma(1,6), est_ers79.conv_exact(1:8), aller79]
est_ers79.conv_cheap
[est_ers79.conv_cheap(1:8)*2*all_gamma(1,6), est_ers79.conv_exact(1:8), aller79(1:8)]
[est_ers99.conv_cheap(1:5)*2*all_gamma(1,8), est_ers99.conv_exact(1:5), aller99(1:5)]
figure
plot(log10(allmv7i2), log10(aller7i2), '-o')
plot(log10(allmv7i2(1:7)), log10(aller7i2(1:7)), '-o')
hold on
xlabel('log(matvecs)')
ylabel('log(error)')
plot(log10(allmv7i2), log10(est_ers7i2.texp_exact), '-o')
plot(log10(allmv7i2), log10(est_ers7i2.fU), '-o')
plot(log10(allmv7i2), log10(est_ers7i2.conv_exact), '-o')
[allNt93, allmv93, aller93, est_ers93] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 3, 200, 8, 1);
[allNt93, allmv93, aller93, est_ers93] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 3, 150, 8, 1);
[allNt93, allmv93, aller93, est_ers93] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 3, 300, 8, 1);
hold on
plot(log10(allmv93), log10(est_ers93.texp), '-o')
plot(log10(allmv93), log10(est_ers93.texp_exact), '-o')
plot(log10(allmv93), log10(est_ers93.fU), '-o')
plot(log10(allmv93), log10(est_ers93.conv_exact), '-o')
polyfit(log10(allmv93(1:5)), log10(aller93(1:5)), 1)
polyfit(log10(allmv93(1:5)), log10(est_ers93.fU(1:5)), 1)
polyfit(log10(allmv93(1:5)), log10(est_ers93.conv_exact(1:5)), 1)
[est_ers93.conv_cheap(1:5)*2*all_gamma(1,8), est_ers93.conv_exact(1:5)]
[est_ers93.conv_cheap(1:5)*2*all_eta(1,11), est_ers93.conv_exact(1:5)]
save errors_harmonic
[allNt99i2, allmv99i2, aller99i2, est_ers99i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 9, 200, 7, 2);
hold on
plot(log10(allmv99i2), log10(est_ers99i2.texp_exact), '-o')
plot(log10(allmv99i2), log10(est_ers99i2.fU), '-o')
plot(log10(allmv99i2), log10(est_ers99i2.conv_exact), '-o')
[allNt9_13i2, allmv9_13i2, aller9_13i2, est_ers9_13i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 13, 150, 7, 2);
[allNt9_13i2, allmv9_13i2, aller9_13i2, est_ers9_13i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 13, 100, 7, 2);
[allNt9_13i2, allmv9_13i2, aller9_13i2, est_ers9_13i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 13, 70, 7, 2);
log10(70)
10^(log10(70) + 0.1)
[allNt9_13i2, allmv9_13i2, aller9_13i2, est_ers9_13i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 13, 90, 7, 2);
[allNt9_13i2, allmv9_13i2, aller9_13i2, est_ers9_13i2] = error_decaySG1(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase(:, end), [0 10], 9, 13, 100, 7, 2);
hold on
plot(log10(allmv9_13i2), log10(est_ers9_13i2.texp_exact), '-o')
hold on
plot(log10(allmv9_13i2), log10(est_ers9_13i2.fU), '-o')
plot(log10(allmv9_13i2), log10(est_ers9_13i2.texp_exact), '-o')
plot(log10(allmv9_13i2), log10(est_ers9_13i2.conv_exact), '-o')
polyfit(log10(allmv9_13i2(1:5)), log10(est_ers9_13i2.conv_exact(1:5)), 1)
polyfit(log10(allmv9_13i2(1:5)), log10(aller9_13i2(1:5)), 1)
[est_ers9_13i2.conv_cheap(1:5)*2*all_eta(1,11), est_ers9_13i2.conv_exact(1:5)]
[est_ers9_13i2.conv_cheap(1:5)*2*all_gamma(1,8), est_ers9_13i2.conv_exact(1:5)]
[est_ers9_13i2.conv_cheap(1:5)*2*all_gamma(2,8), est_ers9_13i2.conv_exact(1:5)]
[est_ers9_13i2.conv_cheap(1:5)*2*all_gamma(2,8), est_ers9_13i2.conv_exact(1:5), aller9_13i2(1:5)]
figure
plot(log10(allmv9_13i2(1:7)), log10(aller9_13i2(1:7)), '-o')
hold on
plot(log10(allmv9_13i2), log10(est_ers9_13i2.conv_exact), '-o')
plot(log10(allmv9_13i2), log10(est_ers9_13i2.conv_cheap), '-o')
plot(log10(allmv9_13i2), log10(est_ers9_13i2.conv_cheap)*2*all_gamma(2,8), '-o')
plot(log10(allmv9_13i2), log10(est_ers9_13i2.conv_cheap*2*all_gamma(2,8)), '-o')
figure
plot(log10(allmv93), log10(est_ers93.conv_exact), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
plot(log10(allmv93), log10(est_ers93.conv_cheap), '-o')
plot(log10(allmv93), log10(est_ers93.conv_exact), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
hold on
plot(log10(allmv93), log10(est_ers93.conv_cheap), '-o')
plot(log10(allmv93), log10(est_ers93.conv_cheap*2*all_eta(1,11)), '-o')
figure
plot(log10(allmv99), log10(aller99), '-o')
hold on
plot(log10(allmv99), log10(est_ers99.conv_exact), '-o')
plot(log10(allmv99), log10(est_ers99.conv_cheap), '-o')
plot(log10(allmv99), log10(est_ers99.conv_cheap*2*all_gamma(1,8)), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
figure
plot(log10(allmv79), log10(aller79), '-o')
hold on
plot(log10(allmv79), log10(est_ers79.conv_exact), '-o')
plot(log10(allmv79), log10(est_ers79.conv_cheap), '-o')
plot(log10(allmv79), log10(est_ers79.conv_cheap*2*all_gamma(1, 6)), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')