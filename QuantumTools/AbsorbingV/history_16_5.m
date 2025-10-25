load('Vabs14_5_18.mat')
[Vopt2c, performance2c, gradmin2c, niter2c, nfevals2c, dif_p2c, dif_V2c, conv2c, alpha_last2c, invHess2c] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 97, 0), Vopt2b, options);
performance2c
[allT2b, allR2b, allnorm2b, allgradnorm2b] = allTRcoef1(Vx2b, [0 40], [0 5], 101);
[allT649, allR649, allnorm649, allgradnorm649] = allTRcoef1(Vopt649x, [0 40], [0 5], 101);
figure
plot(0:0.05:5, allnorm2b)
plot(0:0.05:5, log10(allnorm2b))
hold on
plot(0:0.05:5, log10(allnorm649))
[Vopt2d, performance2d, gradmin2d, niter2d, nfevals2d, dif_p2d, dif_V2d, conv2d, alpha_last2d, invHess2d] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 64, 193, 0), Vopt2c, options);
10^-9.157
ans/3.6886e-10
97/193
193/97
performance2d
figure
Vx2c = Vopt2Vx0lb3b(Vopt2c, 64);
Vx2d = Vopt2Vx0lb3b(Vopt2d, 64);
figure
plot(x, real(Vx2c))
hold on
plot(x, imag(Vx2c))
plot(x, real(Vx2b))
plot(x, imag(Vx2b))
plot(x, real(Vx2d))
plot(x, imag(Vx2d))
[Vopt2e, performance2e, gradmin2e, niter2e, nfevals2e, dif_p2e, dif_V2e, conv2e, alpha_last2e, invHess2e] = quasiNewton(@(Vk) percosV0lb_grad3b(Vk, [0 40], [0.2 5], 128, 97, 0), Vopt2c, options);
Vx_params = Vopt2Vx0lb3b(xmin, 64);
figure
plot(0:40/128:40, Vopt2Vx0lb3b(xmin, 128))
hold on
plot(0:40/128:40, imag(Vopt2Vx0lb3b(xmin, 128)))
percosV0lb_grad3b(xmin, [0 40], [0.2 5], 128, 97, 0)
percosV0lb_grad3b(xmin, [0 40], [0.2 5], 256, 97, 0)
percosV0lb_grad3b(xmin, [0 40], [0.2 5], 512, 97, 0)
percosV0lb_grad3b(xmin, [0 40], [0.2 5], 64, 97, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 64, 193, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 64, 97, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 128, 97, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 128, 193, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 256, 97, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 512, 97, 0)
percosV0lb_grad3b(Vopt2e, [0 40], [0.2 5], 1024, 97, 0)
save Vabs14_5_18