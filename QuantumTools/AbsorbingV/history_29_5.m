[Vabs240b, xabs240b, x240b, p240b, K240b, Nx240b, V0240b] = get_prop_vars(Vf, xdomain240, 40/64, Vx2b);
Vf = @(x) 1 - 1./sqrt(1+x.^2);
[Vabs240b, xabs240b, x240b, p240b, K240b, Nx240b, V0240b] = get_prop_vars(Vf, [-240 240], 40/64, Vx2b);
size(Vabs240b)
figure
plot(-280:0.625:280-0.625, Vabs240b)
plot(-280:0.625:(280-0.625), Vabs240b)
size(-280:0.625:(280-0.625))
40/64
size(-240:0.625:(240-0.625))
plot(-240:0.625:(240-0.625), Vabs240b)
hold on
plot(-240:0.625:(240-0.625), imag(Vabs240b))
whos
dif_p2e
dif_V2e
Vx2e = Vopt2Vx0lb3b(Vopt2e, 64);
figure
plot(x, real(Vx2e))
hold on
plot(x, imag(Vx2e))
[Vabs240e, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-240 240], 40/64, Vx2e);
[uc, mniterc, matvecsc, all_est_erc] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
[uc, mniterc, matvecsc, all_est_erc] = SemiGlobalArnoldi_xp(K240, Vabs240c, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
[ue, mnitere, matvecse, all_est_ere] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
mniterc
mnitere
matvecsc
matvecse
[Vabs480e, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vx2e);
[ue2, mnitere2, matvecse2, all_est_ere2] = SemiGlobalArnoldi_xp(K480, Vabs480e, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
matvecse2
mnitere2
figure
plot(x240(65:705), conj(ue2(449:1089,2)).*ue2(449:1089,2)- conj(ue(65:705,2)).*ue(65:705,2))
plot(x240(65:705), conj(uc2(449:1089,2)).*uc2(449:1089,2)- conj(uc(65:705,2)).*uc(65:705,2))
all_est_ere
all_est_ere2
[ue, mnitere, matvecse, all_est_ere] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
all_est_ere
[ue2, mnitere2, matvecse2, all_est_ere2] = SemiGlobalArnoldi_xp(K480, Vabs480e, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
all_est_ere2
hold on
plot(x240(65:705), conj(ue2(449:1089,2)).*ue2(449:1089,2)- conj(ue(65:705,2)).*ue(65:705,2))
[us, mniters, matvecss, all_est_ers] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) 0, [], fi0480b, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
[us, mniters, matvecss, all_est_ers] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) zeros(768,1), [], fi0480b, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
[us, mniters, matvecss, all_est_ers] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) zeros(768,1), [], fi0240, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false
figure
plot(x240(65:705), conj(fi0480b(449:1089,2)).*fi0480b(449:1089,2)- conj(fi0240(65:705,2)).*fi0240(65:705,2))
plot(x240(65:705), conj(fi0480b(449:1089)).*fi0480b(449:1089)- conj(fi0240(65:705)).*fi0240(65:705))
[us, mniters, matvecss, all_est_ers] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) zeros(768,1), [], fi0240, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
mniters
figure
plot(x240(65:705), conj(us(65:705,2)).*us(65:705,2)- conj(fi0240(65:705,2)).*fi0240(65:705,2))
plot(x240(65:705), conj(us(65:705,2)).*us(65:705,2)- conj(fi0240(65:705)).*fi0240(65:705))
x=0:pi/600:pi;
x=0:40/64:40;
x1=0:pi/600:pi;
f1 = cos(2*x1)-1;
f1 = cos(2*x1)-1 + (cos(2*x1)-1)*1i;
[allT, allR, allnorm] = allTRcoef1(f1, [0 pi], [1 60], 60);
figure
plot(1:60, [allT; allR; allnorm])
plot(1:60, allT)
hold on
plot(1:60, allR)
[allT1, allR1, allnorm1] = allTRcoef1(f1(1:10:601), [0 pi], [1 60], 60);
plot(1:60, allT1)
figure
plot(1:60, allT1-allT)
plot(1:60, allR1-allR)
plot(1:60, log10(abs(allnorm1-allnorm)))
[allT2, allR2, allnorm2] = allTRcoef1(f1(1:40:601), [0 pi], [1 60], 60);
figure
plot(1:60, log10(abs(allnorm2-allnorm)))
[allT3, allR3, allnorm3] = allTRcoef1(f1(1:20:601), [0 pi], [1 60], 60);
plot(1:60, log10(abs(allnorm3-allnorm)))
plot(1:60, (abs(allnorm3-allnorm)))
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 64, 49, 0)
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 64, 97, 0)
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 64, 193, 0)
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 128, 97, 0)
Vx2b128 = Vopt2Vx0lb3b(Vopt2b, 128);
[allT2b, allR2b, allnorm2b] = allTRcoef1(Vx2b, [0 40], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm2b)
[allT2b128, allR2b128, allnorm2b128] = allTRcoef1(Vx2b128, [0 40], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm2b128)
[allT2b128a, allR2b128a, allnorm2b128a] = allTRcoef1([Vx2b128; Vx2b128(end)], [0 40+40/128], [0.2 5], 97);
hold on
plot(0.2:0.05:5, allnorm2b128a)
hold on
plot(0.2:0.05:5, allnorm2b128)
figure
plot(0.2:0.05:5, allnorm2b128./allnorm2b)
plot(0.2:0.05:5, log10(allnorm2b128./allnorm2b))
figure
plot(0.2:0.05:5, log10(allT2b128./allT2b))
hold on
plot(0.2:0.05:5, log10(allTRb128./allTRb))
plot(0.2:0.05:5, log10(allRb128./allR2b))
plot(0.2:0.05:5, log10(allR2b128./allR2b))
figure
plot(x, real(Vx2b))
hold on
plot(0:40/128:40, real(Vx2b128))
plot(0:40/128:40, real(Vx2b128)*sqrt(2))
Vx2b128 = Vopt2Vx0lb3b(Vopt2b*sqrt(2), 128);
plot(0:40/128:40, real(Vx2b128))
plot(x, imag(Vx2b))
plot(0:40/128:40, imag(Vx2b128))
percosV0lb_grad3b(Vopt2b*sqrt(2), [0 40], [0.2 5], 128, 97, 0)
percosV0lb_grad3b(Vopt2b, [0 40], [0.2 5], 64, 97, 0)
[allT2b128, allR2b128, allnorm2b128] = allTRcoef1(Vx2b128, [0 40], [0.2 5], 97);
figure
plot(0.2:0.05:5, log10(allR2b128./allR2b))
plot(0.2:0.05:5, log10(allnorm2b128./allnorm2b))
figure
plot(1:60, (abs(allnorm3-allnorm)))
hold on
clf
plot(1:60, (abs(allnorm3-allnorm))./allnorm)
plot(1:60, log10((abs(allnorm3-allnorm))./allnorm))
figure
plot(1:60, allnorm)
figure
plot(0.2:0.05:5, log10(allnorm2b128./allnorm2b))
plot(0.2:0.05:5, (allnorm2b128./allnorm2b))
percosV0lb_grad3b(Vopt2b*4, [0 40], [0.2 5], 256, 97, 0)
percosV0lb_grad3b(Vopt2b*2, [0 40], [0.2 5], 256, 97, 0)
Vx2b256 = Vopt2Vx0lb3b(Vopt2b*2, 256);
figure
plot(x, real(Vx2b))
hold on
plot(x, imag(Vx2b))
plot(0:40/256:40, real(Vx2b256))
plot(0:40/256:40, imag(Vx2b256))
[allT2b256, allR2b256, allnorm2b256] = allTRcoef1(Vx2b256, [0 40], [0.2 5], 97);
figure
plot(0.2:0.05:5, (allnorm2b256./allnorm2b256))
plot(0.2:0.05:5, (allnorm2b256./allnorm2b128))
percosV0lb_grad3b(Vopt2b*sqrt(8), [0 40], [0.2 5], 512, 97, 0)
Vx2b512 = Vopt2Vx0lb3b(Vopt2b*sqrt(8), 512);
[allT2b512, allR2b512, allnorm2b512] = allTRcoef1(Vx2b512, [0 40], [0.2 5], 97);
figure
plot(0.2:0.05:5, (allnorm2b512./allnorm2b256))
[Vopt9, performance9, gradmin9, niter9, nfevals9, dif_p9, dif_V9, conv9, alpha_last9, invHess9] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 97, 1e-12), Vopt2b, options
pi/40
pi/40*64
[Vopt9, performance9, gradmin9, niter9, nfevals9, dif_p9, dif_V9, conv9, alpha_last9, invHess9] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [pi/40 (pi/40)*64], 64, 64, 1e-12), Vopt2b, options);
[Vopt9, performance9, gradmin9, niter9, nfevals9, dif_p9, dif_V9, conv9, alpha_last9, invHess9] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [pi/40 (pi/40)*64], 64, 64, 1e-12), Vgb, options);
Vx9 = Vopt2Vx0lb3b(Vopt9, 64);
figure
plot(x, real(Vx9))
hold on
plot(x, imag(Vx9))
[allT9, allR9, allnorm9] = allTRcoef1(Vx9, [0 40], [pi/40 (pi/40)*64], 64);
figure
plot(pi/40:pi/40:pi/40*64, allT9)
plot(pi/40:pi/40:pi/40*64, log19(allT9))
plot(pi/40:pi/40:pi/40*64, log10(allT9))
hold on
plot(pi/40:pi/40:pi/40*64, log10(allR9))
percosV0lb_grad3b(Vopt9, [0 40], [pi/40 pi/40*64], 64, 64, 0)
percosV0lb_grad3b(Vopt9, [0 40], [pi/40 pi/40*64], 64, 128, 0)
percosV0lb_grad3b(Vopt9, [0 40], [pi/80 pi/40*64], 64, 128, 0)
percosV0lb_grad3b(Vopt9, [0 40], [pi/40 pi/40*64], 64, 127, 0)
percosV0lb_grad3b(Vopt9*sqrt(2), [0 40], [pi/40 pi/40*64], 128, 127, 0)
Vopt9
[Vopt9a, performance9a, gradmin9a, niter9a, nfevals9a, dif_p9a, dif_V9a, conv9a, alpha_last9a, invHess9a] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [pi/40 (pi/40)*64], 64, 127, 1e-12), Vopt9, options);
performance9a
[Vopt9b, performance9b, gradmin9b, niter9b, nfevals9b, dif_p9b, dif_V9b, conv9b, alpha_last9b, invHess9b] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [pi/40 (pi/40)*64], 128, 127, 1e-12), Vopt9a, options);
[Vopt9b, performance9b, gradmin9b, niter9b, nfevals9b, dif_p9b, dif_V9b, conv9b, alpha_last9b, invHess9b] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [pi/40 (pi/40)*64], 128, 127, 1e-12), Vopt9a*sqrt(2), options);
performance9b
Vx9a = Vopt2Vx0lb3b(Vopt9a, 64);
Vx9b = Vopt2Vx0lb3b(Vopt9b, 128);
figure
plot(x, real(Vx9a))
plot(x, imag(Vx9a))
plot(0:40/128:40, real(Vx9b))
plot(0:40/128:40, imag(Vx9b))
[allT9b, allR9b, allnorm9b] = allTRcoef1(Vx9b, [0 40], [pi/40 (pi/40)*64], 127);
figure
plot(pi/40:pi/80:pi/40*64, log10(allR9b))
hold on
plot(pi/40:pi/80:pi/40*64, log10(allT9b))
plot(pi/40:pi/80:pi/40*64, log10(allnorm9b))
[Vabs240_9b, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-240 240], 40/64, Vx9b);
[Vabs480_9b, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vx9b);
[u9b, mniter9b, matvecs9b, all_est_er9b] = SemiGlobalArnoldi_xp(K240, Vabs240_9b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
mniter9b
mniterb
all_est_er9b
[u9b480, mniter9b480, matvecs9b480, all_est_er9b480] = SemiGlobalArnoldi_xp(K480, Vabs480_9b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480, x480, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
[u9b480, mniter9b480, matvecs9b480, all_est_er9b480] = SemiGlobalArnoldi_xp(K480, Vabs480_9b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
[u9b480, mniter9b480, matvecs9b480, all_est_er9b480] = SemiGlobalArnoldi_xp(K480, Vabs480_9b, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
matvecs9b480
matvecs9b
all_est_er9b480
figure
plot(x240(65:705), conj(u9b480(449:1089,2)).*u9b480(449:1089,2)- conj(u9b(65:705,2)).*u9b(65:705,2))
figure
plot(x240(65:705), conj(ue2(449:1089,2)).*ue2(449:1089,2)- conj(ue(65:705,2)).*ue(65:705,2))
figure
plot(0:pi/960:pi/40*64, fft(fftshift(ue2)))
size(ue2)
plot(0:pi/960:pi/40*64, fft(fftshift(ue2(:,2))))
ue2p = fft(ue2(:,2));
size(ue2p)
size(0:pi/960:(pi/0.625-pi/960))
plot(0:pi/960:(pi/0.625-pi/960), ue2p.*conj(ue2p))
pi/40
figure
plot(x480, u9b40.*conj(u9b480))
plot(x480, u9b480.*conj(u9b480))
uep = fft(ue(65:705,2));
figure
uep = fft(ue(:,2));
plot(0:pi/480:pi/0.625-pi/480, uep.*conj(uep))
plot(0:pi/480:pi/0.625-pi/480, fi0240.*conj(fi0240))
plot(0:pi/480:pi/0.625-pi/480, fft(fi0240).*conj(fft(fi0240)))
plot(-pi/480:2*pi/480:pi/0.625-pi/480, fftshift(fft(fi0240)).*conj(fftshift(fft(fi0240))))
plot(-pi/480:2*pi/480:(pi/0.625-pi/480), fftshift(fft(fi0240)).*conj(fftshift(fft(fi0240))))
uep = fftshift(fft(ue(:,2)));
fi0240p = fftshift(fft(fi0240));
plot(-pi/480:2*pi/480:(pi/0.625-pi/480), fi0240p.*conj(fi0240p))
size(-pi/480:2*pi/480:(pi/0.625-pi/480))
plot(-pi/0.625:2*pi/480:(pi/0.625-pi/480), fi0240p.*conj(fi0240p))
plot(-pi/0.625:2*pi/480:(pi/0.625-pi/480), uep.*conj(uep))
We = psi2wigner(ue(:,2));
viewWigner(ue, [-240 240], 1, );
Wrange = getWrange(ue)
viewWigner(ue, [-240 240], 1, Wrange)
whos
figure
plot(x240, V0240)
plot(x480, V0480)
[u0, mniter0, matvecs0, all_est_er0] = SemiGlobalArnoldi_xp(K240, V0240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
mniter0
mniterb
all_est_er0
figure
plot(x240, u0(:,2).*conj(u0(:,2))
plot(x240, u0(:,2).*conj(u0(:,2)))
plot(x240, ue(:,2).*conj(ue(:,2)))
hold on
plot(x240, u0(:,2).*conj(u0(:,2)))
clf
plot(x240, u0(:,2).*conj(u0(:,2))-ue(:,2).*conj(ue(:,2)))
[u0480, mniter0480, matvecs0480, all_est_er0480] = SemiGlobalArnoldi_xp(K480, V0480, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480, x480, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
[u0480, mniter0480, matvecs0480, all_est_er0480] = SemiGlobalArnoldi_xp(K480, V0480, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 9, 1e-8, 10, 20, false);
mniter0480
all_est_er0480
plot(x240(65:705), conj(u0480(449:1089,2)).*u0480(449:1089,2)- conj(ue2(449:1089,2)).*ue2(449:1089,2))
hold on
plot(x240, u0(:,2).*conj(u0(:,2))-ue(:,2).*conj(ue(:,2)))
a0240 = -x240./(1 + x240.^2).^(3/2);
a0480 = -x480./(1 + x480.^2).^(3/2);
[u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:1e3, 1e4, 7, 7, 1e-5, 10, 20, false);
ma = evmiu(u, a0240);
figure
w = 0:pi/1e3:pi;
maw = dctI(ma);
plot(w, maw)
[u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e3, 1e4, 7, 9, 1e-5, 10, 20, false);
[u2, mniter2, matvecs2, all_est_er2] = SemiGlobalArnoldi_xp(K480, Vabs480, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-5, 10, 20, false);
ma = evmiu(u, a0240);
w = 0:pi/1e3:pi/0.2;
maw = dctI(ma);
plot(w, maw)
plot(w(1:1001), maw(1:1001))
ma2 = evmiu(u2, a0480);
maw2 = dctI(ma2);
hold on
plot(w(1:1001), maw2(1:1001))
figure
plot(w(1:1001), maw2(1:1001)-maw(1:1001))
plot(w(1:1001), abs(maw2(1:1001)-maw(1:1001))./abs(maw2(1:1001)))
[u9b480, mniter9b480, matvecs9b480, all_est_er9b480] = SemiGlobalArnoldi_xp(K480, Vabs480_9b, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
[u9b, mniter9b, matvecs9b, all_est_er9b] = SemiGlobalArnoldi_xp(K240, Vabs240_9b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e], 1e4, 7, 9, 1e-8, 10, 20, false);
[u9b, mniter9b, matvecs9b, all_est_er9b] = SemiGlobalArnoldi_xp(K240, Vabs240_9b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
ma9b = evmiu(u9b, a0240);
ma9b480 = evmiu(u9b480, a0480);
maw9b = dctI(ma9b);
maw9b480 = dctI(ma9b480);
figure
plot(w(1:1001), maw9b(1:1001))
hold on
plot(w(1:1001), maw(1:1001))
plot(w(1:1001), maw9b(1:1001)-maw(1:1001))
clf
plot(w(1:1001), maw9b(1:1001)-maw(1:1001))
plot(w(1:1001), abs(maw9b(1:1001)-maw(1:1001))./maw(1:1001))
plot(w(1:1001), abs(maw9b480(1:1001)-maw9b(1:1001))./abs(maw9b480(1:1001)))
plot(w(1:1001), abs(maw9b480-maw9b./abs(maw9b480))
plot(w(1:1001), abs(maw9b480-maw9b)./abs(maw9b480))
plot(w, abs(maw9b480-maw9b)./abs(maw9b480))
figure
plot(w, abs(maw2-maw)./abs(maw2))
sqrt(0.12)
figure
plot(w, abs(maw2-maw))
hold on
plot(w, abs(maw9b480-maw9b))
plot(w, abs(maw9b480-maw2))
0.2^2/2
figure
plot(x, real(Vx2e))
hold on
plot(x, imag(Vx2e))
2*pi/0.2
[ue, mnitere, matvecse, all_est_ere] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
[ue2, mnitere2, matvecse2, all_est_ere2] = SemiGlobalArnoldi_xp(K480, Vabs480e, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
mae = evmiu(ue, a0240);
mae2 = evmiu(ue2, a0480);
mawe = dctI(mae);
mawe2 = dctI(mae2);
figure
plot(w, abs(mawe2-mawe))
5^2/2
hold on
plot(w, abs(maw2-maw))
figure
plot(w, abs(maw2-maw))
figure
plot(w, abs(mawe2-mawe)./abs(mawe2))
hold on
plot(w, abs(maw2-maw)./abs(maw2))
figure
plot(w, abs(maw2-maw)./abs(maw2))
figure
plot(pi/40:pi/80:pi/40*64, log10(allnorm9b))
plot((pi/40:pi/80:pi/40*64).^2/2, log10(allnorm9b))
figure
plot(w, abs(maw9b480-maw9b)./abs(maw9b480))
[allT9b, allR9b, allnorm9b] = allTRcoef1(Vx2e, [0 40], [0.2 5], 97);
figure
plot((0.2:0.05:5).^2/2, allnorm9b)
[allT9b, allR9b, allnorm9b] = allTRcoef1(Vx9b, [0 40], [pi/40 (pi/40)*64], 127);
[allTe, allRe, allnorme] = allTRcoef1(Vx2e, [0 40], [0.2 5], 97);
performance2e
Vx2e = Vopt2Vx0lb3b(Vopt2e/sqrt(2), 64);
figure
plot(x, real(Vx2e))
hold on
plot(x, imag(Vx2e))
[allTe, allRe, allnorme] = allTRcoef1(Vx2e, [0 40], [0.2 5], 97);
plot((0.2:0.05:5).^2/2, allnorme)
Vx2e128 = Vopt2Vx0lb3b(Vopt2e, 128);
hold on
plot(0:40/128:40, real(Vx2e128))
plot(0:40/128:40, imag(Vx2e128))
[allTe, allRe, allnorme] = allTRcoef1(Vx2e128, [0 40], [0.2 5], 97);
plot((0.2:0.05:5).^2/2, allnorme)
[Vabs240e, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-240 240], 40/64, Vx2e);
[ue, mnitere, matvecse, all_est_ere] = SemiGlobalArnoldi_xp(K240, Vabs240e, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
[ue2, mnitere2, matvecse2, all_est_ere2] = SemiGlobalArnoldi_xp(K480, Vabs480e, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
figure
plot(x240(65:705), conj(ue2(449:1089,2)).*ue2(449:1089,2)- conj(ue(65:705,2)).*ue(65:705,2))
plot(x240(65:705), conj(ue2(449:1089,end)).*ue2(449:1089,end)- conj(ue(65:705,end)).*ue(65:705,end))
mae = evmiu(ue, a0240);
mae2 = evmiu(ue2, a0480);
mawe = dctI(mae);
mawe2 = dctI(mae2);
figure
plot(w, abs(mawe2-mawe)./abs(mawe2))
[Vabs480e, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vx2e);
[ue2, mnitere2, matvecse2, all_est_ere2] = SemiGlobalArnoldi_xp(K480, Vabs480e, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
plot(x240(65:705), conj(ue2(449:1089,end)).*ue2(449:1089,end)- conj(ue(65:705,end)).*ue(65:705,end))
mae2 = evmiu(ue2, a0480);
mawe2 = dctI(mae2);
plot(w, abs(mawe2-mawe)./abs(mawe2))
figure
plot(x, real(Vx2e))
hold on
plot(x, imag(Vx2e))
clf
figure
plot(x480, Vabs480e)
plot(x480, imag(Vabs480e))
plot(x240, imag(Vabs240e))
plot(w, abs(mawe2-mawe)./abs(mawe2))
plot(w, abs(mawe2-mawe))
[allTes, allRes, allnormes] = allTRcoef1([Vx2e128; Vx2e128(end-1:-1:1)], [0 80], [0.2 5], 97);
figure
plot((0.2:0.05:5).^2/2, allnormes)
Vx9b = Vopt2Vx0lb3b(Vopt9b, 128);
figure
plot(0:40/128:40, real(Vx9b))
hold on
plot(0:40/128:40, imag(Vx9b))
size(Vabs240_9b)
clf
plot(x240, real(Vabs240_9b))
Vx9b64 = Vopt2Vx0lb3b(Vopt9b/sqrt(2), 64);
figure
plot(0:40/128:40, real(Vx9b))
hold on
plot(0:40/128:40, imag(Vx9b))
plot(x, real(Vx9b64))
plot(x, imag(Vx9b64))
[Vabs240_9b, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-240 240], 40/64, Vx9b64);
[u9b, mniter9b, matvecs9b, all_est_er9b] = SemiGlobalArnoldi_xp(K240, Vabs240_9b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
mniter9b
all_est_er9b
[Vabs480_9b, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vx9b64);
[u9b480, mniter9b480, matvecs9b480, all_est_er9b480] = SemiGlobalArnoldi_xp(K480, Vabs480_9b, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
figure
plot(x240(65:705), conj(u9b480(449:1089,end)).*u9b480(449:1089,end)- conj(u9b(65:705,end)).*u9b(65:705,end))
figure
plot(0.2:0.05:5, allnorm9b)
plot(pi/40:pi/80:(pi/40)*64, allnorm9b)
percosV0lb_grad3b(Vopt9b, [0 40], [pi/40 pi/40*64], 128, 127, 0)
plot((pi/40:pi/80:(pi/40)*64).^2/2, allnorm9b)
whos
clear We
maw9b = dctI(ma9b);
maw9b480 = dctI(ma9b480);
figure
plot(w, abs(mawe9b480-maw9b))
plot(w, abs(maw9b480-maw9b))
plot(w, abs(maw9b480-maw9b)./abs(maw9b480))
figure
plot(w, abs(mawe2-mawe))
hold on
plot(w, abs(maw9b480-maw9b))
figure
plot(w, abs(maw9b480-maw9b)./abs(maw9b480))
hold on
plot(w, abs(mawe2-mawe)./abs(mawe2))
all_est_er9b480
whos
clear ue ue2 u u2 u9b u9b480
figur
figure
plot((pi/40:pi/80:(pi/40)*64).^2/2, allnorm9b)
plot((pi/40:pi/80:(pi/40)*64), allnorm9b)
plot((pi/40:pi/80:(pi/40)*64).^2/2, allnorm9b)
figure
plot(w, abs(maw9b480-maw9b)./abs(maw9b480))
figure
plot(w, maw9b)
plot(w, maw)
plot(w, maw9b)
plot(w(1:1001), maw9b(1:1001))
hold on
plot(w(1:1001), maw9b480(1:1001))
[u0, mniter0, matvecs0, all_est_er0] = SemiGlobalArnoldi_xp(K240, V0240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
[u0480, mniter0480, matvecs0480, all_est_er0480] = SemiGlobalArnoldi_xp(K480, V0480, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, 0:0.2:1e3, 1e4, 7, 9, 1e-8, 10, 20, false);
ma0 = evmiu(u0, a0240);
ma0480 = evmiu(u0480, a0480);
ma0w = dctI(ma0);
ma0480w = dctI(ma0480);
figure
plot(w, ma0480)
plot(w, ma0480w)
hold on
plot(w, mawe)
clf
plot(w(1:1001), ma0480w(1:1001))
hold on
plot(w(1:1001), ma0w(1:1001))
plot(w(1:1001), ma0480w(1:1001))
clear u0 u0480
whos
save Vabs29_5_18