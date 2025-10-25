load coulomb_optV240
whos
Vf = @(x)1-1./sqrt(x.^2+1)
[~, ~, ~, ~, P240] = gsV(Vf, [-240 240], Nx240);
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
fieldwg = 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015);
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(fieldwg)/dctfactor;
figure
plot(t, fieldtg)
figure
plot(w, fieldwg)
plot(w(1:101), fieldwg(1:101))
[~, ~, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
(0.05/5.33e-9)^2
5.33e-9*sqrt(1e14)
5.33e-9*sqrt(1e15)
figure
plot(w(1:1001), evawg(1:1001))
npsig=sqnorm(psig);
figure
plot(t, npsig)
[~, ~, psige, evatge, evawge, evmiutge, evmiuwge, mnitercge, Jge, J1ge, J2ge, Jorthge, Jpnormge] = guessresults_pnaE0b(sqrt(0.95)*P240(:,1) + sqrt(0.05)*P240(:,2), Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jge, J1ge, J2ge, Jorthge, Jpnormge]
npsige=sqnorm(psige);
figure
plot(t, npsige)
figure
plot(w(1:1001), evawge(1:1001))
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
save field78poster fieldt fieldw evat evaw evmiut evmiuw relE conv1 niter1 mallniterc J1 maxgrad weight
figure
plot(w(1:1001), evaw(1:1001))
figure
plot(t, fieldt)
save field78poster fieldt fieldw evat evaw evmiut evmiuw relE conv1 niter1 mallniterc J1 maxgrad weight fieldtg fieldwg evatg evawg evmiutg evmiuwg mnitercg Jg J1g J2g Jorthg Jpnormg
hold on
plot(t, fieldtg, 'r')
figure
plot(w(1:1001), evaw(1:1001))
%-- 22/04/2015 12:48 --%
load('field78poster.mat')
load coulomb_optV240
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(t, fieldt)
dctfactor = 1e3/(sqrt(5e3*pi))
figure
[~, ~, psi] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
npsi=sqnorm(psi);
figure
plot(t, npsi)
figure
plot(w(1:101), fieldw(1:101))
fieldwg = 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015);
fieldtg = dctI(fieldwg)/dctfactor;
hold on
plot(w(1:101), fieldwg(1:101), 'r')
iwfield = instwcos(fieldt, 1e3);
figure
plot(t, iwfield)
figure
plot(w(1:1001), evaw(1:1001))
300/6
Epong = 0.05^2/(4*0.06)
3.17*Epong + 1
[fieldth, fieldwh, psih, evath, evawh, evmiuth, evmiuwh, relEh, convh, niterh, mallniterch, J1h, maxgradh, weighth] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 1.5, 3), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(w(1:1001), evawh(1:1001))
npsih=sqnorm(psih);
figure
plot(t, npsih)
figure
plot(t, fieldth)
208/0.06
2.08/0.06
2.06/0.06
2.1/0.06
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35, niter35, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 2.09, 2.11), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35, niter35, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 2.09, 2.11), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35, niter35, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 2.09, 2.11), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35, niter35, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh, @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 2.09, 2.11), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35, niter35, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh, @(w) 1e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 2.09, 2.11), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35, niter35, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh, @(w) 1e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 2.09, 2.11), 0, 1e-2, 1e3, 0.2, 9, 9, 1e-3, 1e3);
save field1530 fieldth fieldwh evath evawh evmiuth evmiuwh relEh convh niterh mallniterch J1h maxgradh weighth
%%%%%%%%%%%%%%
[fieldts, fieldws, psis, evats, evaws, evmiuts, evmiuws, relEs, convs, niters, mallnitercs, J1s, maxgrads, weights] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*(exp(-(w-0.06).^2/(2*0.01^2))+exp(-(w-0.03).^2/(2*0.005^2))), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtd, fieldwd, psid, evatd, evawd, evmiutd, evmiuwd, relEd, convd, niterd, mallnitercd, J1d, maxgradd, weightd] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*(exp(-(w-0.06).^2/(2*0.01^2))+exp(-(w-0.12).^2/(2*0.02^2))), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);