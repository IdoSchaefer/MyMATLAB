whos
clear all
whos
Vabs_efficiency1(Vopt649x, [0 40], [0.2 5], Nk)
Vabs_efficiency1(Vopt649x, [0 40], [0.2 5], 49)
0.625*64
pi/0.625
pi/40
pi/240
Vabs_efficiency1(Vopt649x, [0 40], [0.1 5], 50)
2*pi/40
Vabs_efficiency1(Vopt649x, [0 64], [0.2 5]/0.625, 49)
Vabs_efficiency1(Vopt649x, [0 64], [0.2 5]*0.625, 49)
pi/1
pi/0.625
5*0.625
size(Vopt649x)
Vabs_efficiency1(Vopt649x, [0 64], [0.2*40/64 5*0.625], 49)
40/64
TRcoef1(0.2, Vopt649x, [0 64], 65)
TRcoef1(0.2*0.625, Vopt649x, [0 64], 65)
TRcoef1(5*0.625, Vopt649x, [0 64], 65)
TRcoef1(3*0.625, Vopt649x, [0 64], 65)
[T, R] = TRcoef1(0.2*0.625, Vopt649x, [0 64], 65)
[T, R] = TRcoef1(0.2, Vopt649x, [0 40], 65)
[T, R] = TRcoef1(0.2, real(Vopt649x), [0 40], 65)
[T, R] = TRcoef1(0.2*0.625, real(Vopt649x), [0 64], 65)
L = 16*sqrt(pi);
Nx = 128;
dx = L/Nx;
x = (-L/2:dx:(L/2 - dx)).';
% Constructing the kinetic energy matrix diagonal in the p domain:
p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
% The potential energy matrix diagonal in the x domain:
V = x.^2/2;
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
[U, mniter, matvecs] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*1e-3*sin(pi*t/100).^2.*cos(5*t), v), @(u1, t1, u2, t2) -1i*x*(1e-3*sin(pi*t1/100).^2.*cos(5*t1) - ones(1, Nt_ts)*1e-3*sin(pi*t2/100).^2.*cos(5*t2)).*u1, 0, [], [-188*1i, 1i], fi0, 0:0.1:100, 1e3, 7, 7, 1e-5);
[U, mniter, matvecs] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*1e-3*sin(pi*t/100).^2.*cos(5*t), v), @(u1, t1, u2, t2) -1i*x*(1e-3*sin(pi*t1/100).^2.*cos(5*t1) - ones(1, 7)*1e-3*sin(pi*t2/100).^2.*cos(5*t2)).*u1, 0, [], [-188*1i, 1i], fi0, 0:0.1:100, 1e3, 7, 7, 1e-5);
[U, mniter, matvecs] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*1e-3*sin(pi*t/100).^2.*cos(5*t), v), @(u1, t1, u2, t2) -1i*x*(1e-3*sin(pi*t1/100).^2.*cos(5*t1) - ones(1, 7)*1e-3*sin(pi*t2/100).^2.*cos(5*t2)).*u1, 0, [], [-188*1i, 1i], fi0, 0:0.1:100, 2e3, 7, 7, 1e-5);
mniter
[U, mniter, matvecs] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*1e-3*sin(pi*t/100).^2.*cos(5*t), v), @(u1, t1, u2, t2) -1i*x*(1e-3*sin(pi*t1/100).^2.*cos(5*t1) - ones(1, 7)*1e-3*sin(pi*t2/100).^2.*cos(5*t2)).*u1, 0, [], [-188*1i, 1i], fi0, 0:0.1:100, 2e3, 7, 9, 1e-5);
mniter
norm(U(:,end))
matvecs
[U1, mniter1, matvecs1] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*1e-3*sin(pi*t/100).^2.*cos(5*t), v), @(u1, t1, u2, t2) -1i*x*(1e-3*sin(pi*t1/100).^2.*cos(5*t1) - ones(1, 7)*1e-3*sin(pi*t2/100).^2.*cos(5*t2)).*u1, 0, [], [-188*1i, 1i], fi0, 0:0.1:100, 2e4, 9, 13, eps);
[U1, mniter1, matvecs1] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*1e-3*sin(pi*t/100).^2.*cos(5*t), v), @(u1, t1, u2, t2) -1i*x*(1e-3*sin(pi*t1/100).^2.*cos(5*t1) - ones(1, 9)*1e-3*sin(pi*t2/100).^2.*cos(5*t2)).*u1, 0, [], [-188*1i, 1i], fi0, 0:0.1:100, 2e4, 9, 13, eps);
mniter1
matvecs1
norm(U(:,end) - U1(:,end))
viewVPmiux(U, @(x) x.^2, x, 1e-3*sin(pi*t/100).^2.*cos(5*t), x, 0.1)
t=0:0.1:100;
viewVPmiux(U, @(x) x.^2, x, 1e-3*sin(pi*t/100).^2.*cos(5*t), x, 0.1)
figure
plot(t, fi0.'*U)
plot(t, (fi0.'*U).*conj(fi0.'*U))
plot(t, 1-(fi0.'*U).*conj(fi0.'*U))
2*pi/5
ans/2
options = default_op_qn
options.minimal_f = 0;
options.plot = @plot;
Vg = rand(32, 1)-0.5
pi/0.625
5-0.2
[Vopt, performance, gradmin, niter, nfevals, dif_p, dif_V, conv, alpha_last, invHess] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg, options);
Vx = Vopt2Vx0lb3(Vopt, 64);
x=0:40/64:64;
plot(x, real(Vx)
figure
plot(x, real(Vx))
size(Vx)
x=0:40/64:40;
plot(x, real(Vx))
hold on
plot(x, imag(Vx))
options_fminunc = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', 1e-8, 'TolX', 1e-8, 'MaxIter', 1e4, 'MaxFunEvals', 1e5);
[Vopt1, optval1, flag1, data1] = fminunc(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg, options_fminunc);
Vg1 = rand(32, 1)-0.5;
[Vopt1, optval1, flag1, data1] = fminunc(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg1, options_fminunc);
options_fminunc = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', 1e-14, 'TolX', 1e-14, 'MaxIter', 1e4, 'MaxFunEvals', 1e5);
[Vopt1, optval1, flag1, data1] = fminunc(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vopt1, options_fminunc);
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg1, options);
options.plot = @(x,y) plot(x, log10(y));
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg1, options);
Vx1 = Vopt2Vx0lb3(Vopt1, 64);
figure
plot(x, real(Vx1))
hold on
plot(x, imag(Vx1))
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg1, options);
175/80
650/330
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg1, options);
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-10), Vg1, options);
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-9), Vg1, options);
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-8), Vg1, options);
Vopt1b
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-7), Vg1, options);
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-6), Vg1, options);
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-5), Vg1, options);
Vopt1b
[Vopt1b, performance1, gradmin1, niter1, nfevals1, dif_p1, dif_V1, conv1, alpha_last1, invHess1] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 1e-4), Vg1, options);
[Vopt1b_2, performance1_2, gradmin1_2, niter1_2, nfevals1_2, dif_p1_2, dif_V1_2, conv1_2, alpha_last1_2, invHess1_2] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vg1b, options);
[Vopt1b_2, performance1_2, gradmin1_2, niter1_2, nfevals1_2, dif_p1_2, dif_V1_2, conv1_2, alpha_last1_2, invHess1_2] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vopt1b, options);
performance1b_2
performance1_2
performance
options1_2=options
options1.f_termination = @(dif_f,dif_x,fmin,gradmin,xmin)terminate_reldif(dif_f,dif_x,fmin,gradmin,xmin,1e-14,1e-14)
clear options1
options1_2.f_termination = @(dif_f,dif_x,fmin,gradmin,xmin)terminate_reldif(dif_f,dif_x,fmin,gradmin,xmin,1e-14,1e-14)
options1_2.invHess0 = invHess1_2
[Vopt1b_3, performance1_3, gradmin1_3, niter1_3, nfevals1_3, dif_p1_3, dif_V1_3, conv1_3, alpha_last1_3, invHess1_3] = quasiNewton(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vopt1b_2, options1_2);
Vopt1b_2
Vx1b = Vopt2Vx0lb3(Vopt1b, 64);
Vx1b_2 = Vopt2Vx0lb3(Vopt1b_2, 64);
figure
plot(x, real(Vx1b))
hold on
plot(x, imag(Vx1b))
figure
plot(x, real(Vx1b_2))
hold on
plot(x, imag(Vx1b_2))
[Vopt1b_3, optval1b_3, flag1b_3, data1b_3] = fminunc(@(Vk) percosV0lb_grad3(Vk, [0 40], [0.2 5], 64, 49, 0), Vopt1b, options_fminunc);
performance1_2
performance1_2-optval1b_3
Vopt1b_2-Vopt1b_3
[Vopt2, performance2, gradmin2, niter2, nfevals2, dif_p2, dif_V2, conv2, alpha_last2, invHess2] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 49, 0), Vg, options);
performance2
performance1
performance
performance1_2
performance1_2-performance2
niter
niter2
Vx2 = Vopt2Vx0lb3b(Vopt2, 64);
Vx2-Vx1_2
Vx2-Vx1b_2
Vopt1b_2-Vopt2
Vgb = [Vg(2:16); -(Vg(1)/2 + sum(Vg(2:16))); Vg(18:32); -(Vg(17)/2 + sum(Vg(18:32)))];
Vxg = Vopt2Vx0lb3b(Vg, 64);
Vxg = Vopt2Vx0lb3(Vg, 64);
Vxgb = Vopt2Vx0lb3b(Vgb, 64);
Vxg-Vxgb
[Vopt2b, performance2b, gradmin2b, niter2b, nfevals2b, dif_p2b, dif_V2b, conv2b, alpha_last2b, invHess2b] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 49, 0), Vgb, options);
Vopt2b
performance2b
performance2
niter
niter2
niter2b
Vx2b = Vopt2Vx0lb3b(Vopt2b, 64);
figure
plot(x, real(Vx2b))
hold on
plot(x, imag(Vx2b))
Vopt2b
Vopt2
Vopt1_2
Voptb1_2
Vopt1b_2
DVx2b = Dcosines(Vx2b, 40);
D2Vx2b = D2cosines(Vx2b, 40);
figure
plot(x, real(DVx2b))
plot(x, imag(DVx2b))
plot(x, real(D2Vx2b))
plot(x, imag(D2Vx2b))
D3Vx2b = Dcosines(D2Vx2b, 40);
plot(x, real(D3Vx2b))
plot(x, imag(D3Vx2b))
Vg5 = [Vgb(1:16); abs(Vgb(17:32))]
Vg5 = rand(32, 1)-0.5;
Vg5(17:32) = sqrt(abs(Vg5(17:32)))
[Vopt3, performance3, gradmin3, niter3, nfevals3, dif_p3, dif_V3, conv3, alpha_last3, invHess3] = quasiNewton(@(Vk) percosV0lb_grad5(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5, options);
Vopt3
figure
Vx3 = Vopt2Vx0lb5(Vopt3, 64);
plot(x, real(Vx3))
hold on
plot(x, imag(Vx3))
Vg5_1 = rand(32, 1)-0.5;
Vg5_1(17:32) = sqrt(abs(Vg5_1(17:32)));
[Vopt3_1, performance3_1, gradmin3_1, niter3_1, nfevals3_1, dif_p3_1, dif_V3_1, conv3_1, alpha_last3_1, invHess3_1] = quasiNewton(@(Vk) percosV0lb_grad5(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5_1, options);
Vopt3_1
Vg5_3 = rand(32, 1)-0.5;
Vg5_2 = rand(32, 1)-0.5;
clear Vg5_3
[Vopt3_2, performance3_2, gradmin3_2, niter3_2, nfevals3_2, dif_p3_2, dif_V3_2, conv3_2, alpha_last3_2, invHess3_2] = quasiNewton(@(Vk) percosV0lb_grad5(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5_2, options);
Vopt3_2
[Vopt3m, optval3m, flag3m, data3m] = fminunc(@(Vk) percosV0lb_grad5(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5, options_fminunc);
performance3
[Vopt4, performance4, gradmin4, niter4, nfevals4, dif_p4, dif_V4, conv4, alpha_last4, invHess4] = quasiNewton(@(Vk) percosV0lb_grad6(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5, options);
performance4
figure
Vopt4
Vx4 = Vopt2Vx0lb6(Vopt4, 64);
[Vopt4m, optval4m, flag4m, data4m] = fminunc(@(Vk) percosV0lb_grad6(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5, options_fminunc);
performance4
figure
plot(x, real(Vx4))
hold on
plot(x, imag(Vx4))
performance2b
Vx2b.*conj(Vx2b)
sum(Vx2b.*conj(Vx2b))
[Vopt2bp, performance2bp, gradmin2bp, niter2bp, nfevals2bp, dif_p2bp, dif_V2bp, conv2bp, alpha_last2bp, invHess2bp] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 49, 1e-12), Vopt2b, options);
figure
percosV0lb_grad3b(Vopt3bp, [0 40], [0.2 5], 64, 49, 1e-12)
percosV0lb_grad3b(Vopt2bp, [0 40], [0.2 5], 64, 49, 1e-12)
percosV0lb_grad3b(Vopt2bp, [0 40], [0.2 5], 64, 49, 0)
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 64, 49, 0)
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 64, 49, 1e-12)
figure
Vx2bp = Vopt2Vx0lb3b(Vopt2bp, 64);
plot(x, real(Vx2bp))
hold on
plot(x, imag(Vx2bp))
hold on
plot(x, real(Vx2b))
plot(x, imag(Vx2b))
[Vopt2bp1, performance2bp1, gradmin2bp1, niter2bp1, nfevals2bp1, dif_p2bp1, dif_V2bp1, conv2bp1, alpha_last2bp1, invHess2bp1] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 49, 1e-12), Vgb, options);
performance2bp1 - performance2bp
nfevals4
niter4
data4m
[Vopt5, performance5, gradmin5, niter5, nfevals5, dif_p5, dif_V5, conv5, alpha_last5, invHess5] = quasiNewton(@(Vk) percosV0lb_grad7(Vk, [0 40], [0.2 5], 64, 49, 0), Vg5, options);
[Vopt5, performance5, gradmin5, niter5, nfevals5, dif_p5, dif_V5, conv5, alpha_last5, invHess5] = quasiNewton(@(Vk) percosV0lb_grad7(Vk, [0 40], [0.2 5], 64, 49, 0), Vg, options);
figure
plot(0:0.625:40, Vximag)
plot(0:0.625:40, Vxreal)
((2:2:2*Nk_free).').*Vk((Nk_free + 1):2:3*Nk_free)
[Vopt5, performance5, gradmin5, niter5, nfevals5, dif_p5, dif_V5, conv5, alpha_last5, invHess5] = quasiNewton(@(Vk) percosV0lb_grad7(Vk, [0 40], [0.2 5], 64, 49, 0), Vg, options);
((2:2:2*Nk_free).').*Vk((Nk_free + 1):2:3*Nk_free)
((2:2:2*Nk_free).').*Vk((Nk_free + 2):2:3*Nk_free)
performance5
Vgbar = rand(64,1)-0.5;
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(V) Vabs_efficiency1(V, [40/64 40], [0.2 5], 49), Vgbar, options);
Vgbar = rand(128,1)-0.5;
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(V) Vabs_efficiency1(V, [40/64 40], [0.2 5], 49), Vgbar, options);
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1(Vk, [0 40], [0.2 5], 64, 49k, 0, 0, 0), Vg, options);
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg, options);
Vg
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg1, options);
whos
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg, options);
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg1, options);
Vg6 = [rand(16, 1)-0.5;rand(16,1)];
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1b(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg1, options);
[Vopt6, performance6, gradmin6, niter6, nfevals6, dif_p6, dif_V6, conv6, alpha_last6, invHess6] = quasiNewton(@(Vk) percosV0lb_grad1b(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg6, options);
performance6
figure
Vx6 = Vopt2Vx0lb1b(Vopt6, 64);
plot(x, real(Vx6))
hold on
plot(x, imag(Vx6))
[Vopt6_1, performance6_1, gradmin6_1, niter6_1, nfevals6_1, dif_p6_1, dif_V6_1, conv6_1, alpha_last6_1, invHess6_1] = quasiNewton(@(Vk) percosV0lb_grad1b(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg6_1, options);
Vg6_1 = [rand(16, 1)-0.5;rand(16,1)];
[Vopt6_1, performance6_1, gradmin6_1, niter6_1, nfevals6_1, dif_p6_1, dif_V6_1, conv6_1, alpha_last6_1, invHess6_1] = quasiNewton(@(Vk) percosV0lb_grad1b(Vk, [0 40], [0.2 5], 64, 49, 0, 0, 0), Vg6_1, options);
performance6_1
figure
Vx6_1 = Vopt2Vx0lb1b(Vopt6_1, 64);
plot(x, real(Vx6_1))
hold on
plot(x, imag(Vx6_1))
[Vopt7, performance7, gradmin7, niter7, nfevals7, dif_p7, dif_V7, conv7, alpha_last7, invHess7] = quasiNewton(@(V) perVabs(V, [0 40], [0.2 5], 0), Vg7, options);
Vg7 = [rand(64,1)-0.5; -rand(64,1)];
[Vopt7, performance7, gradmin7, niter7, nfevals7, dif_p7, dif_V7, conv7, alpha_last7, invHess7] = quasiNewton(@(V) perVabs(V, [0 40], [0.2 5], 0), Vg7, options);
[Vopt7, performance7, gradmin7, niter7, nfevals7, dif_p7, dif_V7, conv7, alpha_last7, invHess7] = quasiNewton(@(V) perVabs(V, [0 40], [0.2 5], 49, 0), Vg7, options);
figure
plot(1:64, xmin(1:64))
hold on
plot(1:64, xmin(65:128))
performance7
figure
plot(x, real(Vopt7))
plot(0:40/63:40, real(Vopt7))
size(Vopt7)
plot(0:40/63:40, Vopt7(1:64))
hold on
plot(0:40/63:40, Vopt7(65:128))
[Vopt8, performance8, gradmin8, niter8, nfevals8, dif_p8, dif_V8, conv8, alpha_last8, invHess8] = quasiNewton(@(V) perVabs(V, [40/16 40], [0.2 5], 49, 0), Vg8, options);
Vg1
[Vopt8, performance8, gradmin8, niter8, nfevals8, dif_p8, dif_V8, conv8, alpha_last8, invHess8] = quasiNewton(@(V) perVabs(V, [40/16 40], [0.2 5], 49, 0), Vg1, options);
Vg8 = [rand(16,1)-0.5; -rand(16,1)];
[Vopt8, performance8, gradmin8, niter8, nfevals8, dif_p8, dif_V8, conv8, alpha_last8, invHess8] = quasiNewton(@(V) perVabs(V, [40/16 40], [0.2 5], 49, 0), Vg8, options);
figure
plot(40/16:40/16:40, Vopt8(1:16))
hold on
plot(40/16:40/16:40, Vopt8(17:32))
[T, R, total_norm] = TRcoef1(0.1, -100*i, [0,1], 1)
[T, R, total_norm] = TRcoef1(0.1, -100*1i, [0,1], 1)
[T, R, total_norm] = TRcoef1(0.1, [-100*1i, -100*1i], [0,1], 2)
[T, R, total_norm] = TRcoef1(0.1, [-1e3*1i, -1e3*1i], [0,1], 2)
[T, R, total_norm] = TRcoef1(0.1, [-1e3*1i, -1e3*1i], [0,100], 2)
[T, R, total_norm] = TRcoef1(0.1, [-1e3*1i, -1e3*1i], [0,10], 2)
[T, R, total_norm] = TRcoef1(0.1, [1e3*1i, 1e3*1i], [0,10], 2)
[T, R, total_norm] = TRcoef1(0.1, [1e2*1i, 1e2*1i], [0,10], 2)
whos
clear U U1
performance2
figure
plot(x, real(Vx2))
hold on
plot(x, imag(Vx2))
performance2b
figure
plot(x, imag(Vx2b))
plot(x, real(Vx2b))
hold on
plot(x, imag(Vx2b))