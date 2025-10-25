test_harmonic
[U50, mniter50, matvecs50, max_errors50, history50] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - ones(1, Nt_ts)*cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, 0:0.2:10, 50, Nt_ts, 17, tol);
fztest
[U50a, mniter50a, matvecs50a, max_errors50a, history50a] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - ones(1, Nt_ts)*cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, 0:0.2:10, 50, Nt_ts, 30, tol);
U50Ea = P\U50a;
abs(f_fun(-1i*195, [0.05:0.05:0.2], 9, 1e-5))*195^9
mx50 = evmiu(U50, x);
mx50a = evmiu(U50a, x);
t50 = 0:0.2:10;
error50 = mx50 - (-0.5*sin(t50).*t50)
error50a = mx50a - (-0.5*sin(t50).*t50);
error50a
figure
plot(t50, conj(U50Ea).*U50Ea)
figure
plot(t50, log10(conj(U50Ea).*U50Ea))
hold on
figure
plot(t50, log10(conj(U50E).*U50E))
U50Ea = P\U50a;
U50E = P\U50;
plot(t50, log10(conj(U50E).*U50E))
plot(t50, log10(conj(U50E(:,1:11)).*U50E(:, 1:11)))
plot(t50(1:11), log10(conj(U50E(:,1:11)).*U50E(:, 1:11)))
figure
plot(t50(1:11), log10(conj(U50Ea(:,1:11)).*U50Ea(:, 1:11)))
element_erU50 = ((conj(U50E).*U50E)-(conj(U50Ea).*U50Ea))./(conj(U50Ea).*U50Ea);
element_erU50(end, 2:end)./element_erU50(end, 1:end-1)
history50.niter
history50a.niter
history50.f_error
history50.fUerror
history50a.fUerror
error50a
for j=0:8, abs(f_fun(-1i*195, [0.05:0.05:0.2], 9, 1e-5)*195^9)./abs(f_fun(-1i*195, [0.05:0.05:0.2], j, 1e-5)*195^j), end
[U, mniter, matvecs, max_errors, history] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - ones(1, Nt_ts)*cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, t, Nts, Nt_ts, Ncheb, tol);
max_errors
max_errors50a
figure
UE = P\U;
plot(t, log10(conj(U).*U))
plot(t, log10(conj(UE).*UE))