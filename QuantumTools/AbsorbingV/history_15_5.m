load('Vabs14_5_18.mat')
load coulomb_optV240
whos
Vf
Vf = @(x) 1-1./sqrt(x.^2 + 1)
figure
plot(x240, Vf(x240))
[Vabs240b, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-240 240], 40/64, Vx2b);
hold on
plot(x240, Vabs240b)
hold on
plot(x240, imag(Vabs240b))
[u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], 5e3, 7, 7, 1e-5, 10, 20, false);
[u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 5e3, 7, 7, 1e-5, 10, 20, false);
[ub, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 5e3, 7, 7, 1e-5, 10, 20, false);
[ub, mniterb, matvecsb, all_est_erb] = SemiGlobalArnoldi_xp(K240, Vabs240b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 5e3, 7, 7, 1e-5, 10, 20, false);
mniterb
matvecs
matvecsb
figure
plot(x240, conj(u(:,2)).*u(:,2))
norm(u(:,2))
norm(u(:,2)).^2
norm(ub(:,2)).^2
performance2b
figure
plot(x240, conj(ub(:,2)).*ub(:,2))
figure
plot(x240, conj(ub(:,2)).*ub(:,2)- conj(u(:,2)).*u(:,2))
x240(65)
x240(384)
x240(385)
385+384
x240(769)
x240(704)
x240(705)
norm(ub(65:705,2)).^2
norm(u(65:705,2)).^2
norm(ub(65:705,2)).^2-norm(u(65:705,2)).^2
figure
plot(x240, conj(ub(65:705,2)).*ub(65:705,2)- conj(u(65:705,2)).*u(65:705,2))
plot(x240(65:705), conj(ub(65:705,2)).*ub(65:705,2)- conj(u(65:705,2)).*u(65:705,2))
plot(x240(65:705), (conj(ub(65:705,2)).*ub(65:705,2)- conj(u(65:705,2)).*u(65:705,2))./(conj(ub(65:705,2)).*ub(65:705,2)))
[Vabs480b, xabs480, x480, p480, K480, Nx480, V0480] = get_prop_vars(Vf, [-480 480], 40/64, Vx2b);
[ub2, mniterb2, matvecsb2, all_est_erb2] = SemiGlobalArnoldi_xp(K480, Vabs480b, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480, x480, [0 1e3], 5e3, 7, 7, 1e-5, 10, 20, false);
fi0480b = gsV(Vf, xdomain480, Nx480);
fi0480b = gsV(Vf, [-480 480], Nx480);
[ub2, mniterb2, matvecsb2, all_est_erb2] = SemiGlobalArnoldi_xp(K480, Vabs480b, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 5e3, 7, 7, 1e-5, 10, 20, false);
[ub2, mniterb2, matvecsb2, all_est_erb2] = SemiGlobalArnoldi_xp(K480, Vabs480b, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
x240(385)
x480(385)
x480(385+768)
385+768
x480(385+704)
x480(429)
x480(430)
x480(449)
385+704
figure
plot(x240(65:705), conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(u(65:705,2)).*u(65:705,2))
hold on
plot(x240(65:705), conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(u(65:705,2)).*u(65:705,2))
mniterb2
[ub, mniterb, matvecsb, all_est_erb] = SemiGlobalArnoldi_xp(K240, Vabs240b, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
[u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
matvecsb
matvecs
figure
plot(x240(65:705), conj(ub(65:705,2)).*ub(65:705,2)- conj(u(65:705,2)).*u(65:705,2))
[Vabs480, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vopt649);
[Vabs480, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vopt649x);
[u2, mniter2, matvecs2, all_est_er2] = SemiGlobalArnoldi_xp(K480, Vabs480, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
figure
plot(x240(65:705), conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(ub(65:705,2)).*u(65:705,2))
plot(x240(65:705), conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(ub(65:705,2)).*ub(65:705,2))
hold on
plot(x240(65:705), conj(u2(449:1089,2)).*u2(449:1089,2)- conj(u(65:705,2)).*u(65:705,2))
figure
plot(x240(65:705), conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(u(65:705,2)).*u(65:705,2))
hold on
plot(x240(65:705), conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(ub(65:705,2)).*ub(65:705,2))
performance2
performance2b
[Vabs240c, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-240 240], 40/64, Vx2);
[uc, mniterc, matvecsc, all_est_erc] = SemiGlobalArnoldi_xp(K240, Vabs240c, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
matvecs
matvecsb
matvecsc
figure
plot(x240(65:705), conj(ub(65:705,2)).*ub(65:705,2)- conj(uc(65:705,2)).*uc(65:705,2))
hold on
plot(x240(65:705), conj(ub(65:705,2)).*ub(65:705,2)- conj(u(65:705,2)).*u(65:705,2))
norm(conj(ub(65:705,2)).*ub(65:705,2)- conj(u(65:705,2)).*u(65:705,2))
norm(conj(ub(65:705,2)).*ub(65:705,2)- conj(uc(65:705,2)).*uc(65:705,2))
norm(conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(u2(65:705,2)).*u2(65:705,2))
norm(conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(u2(449:1089,2)).*u2(449:1089,2))
[uc2, mniterc2, matvecsc2, all_est_erc2] = SemiGlobalArnoldi_xp(K480, Vabs480c, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
[Vabs480c, ~, ~, ~, ~, ~, ~] = get_prop_vars(Vf, [-480 480], 40/64, Vx2);
[uc2, mniterc2, matvecsc2, all_est_erc2] = SemiGlobalArnoldi_xp(K480, Vabs480c, @(u,x,t) -xabs480*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0480b, x480, [0 1e3], 1e4, 7, 7, 1e-5, 10, 20, false);
norm(conj(ub2(449:1089,2)).*ub2(449:1089,2)- conj(uc2(449:1089,2)).*uc2(449:1089,2))
whos
save Vabs14_5_18