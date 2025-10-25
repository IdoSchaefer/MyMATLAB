[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, relE78, conv78, niter78, mallniterc78, J178, maxgrad78, weight78] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1e-2*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1e-2*50*sech(50*(x-0.9))^2, @(w) 20*(exp(-(w-0.06).^2/(2*0.005^2))).*sin((w-0.06)*pi/0.015), @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, relE78, conv78, niter78, mallniterc78, J178, maxgrad78, weight78] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1e-2*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1e-2*50*sech(50*(x-0.9))^2, @(w) 20*(exp(-(w-0.06).^2/(2*0.005^2))).*sin((w-0.06)*pi/0.015), @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldtg)
hold on
plot(t, fieldt78)
plot(t, fieldt90b)
norm(psi78(:,end))
norm(psig(:,end))
figure
plot(w90(1:1001), log10(abs(evawg(1:1001).^2)))
hold on
plot(w90(1:1001), log10(abs(evaw78(1:1001).^2)))
figure
plot(w90(1:51), fieldwg(1:51))
hold on
plot(w90(1:51), fieldw90(1:51))
plot(w90(1:51), fieldw78(1:51))
figure
plot(w90(1:1001), evaw78(1:1001))
plot(w90(1:1001), evawg(1:1001))
hold on
plot(w90(1:1001), evaw78(1:1001))
[fieldt78_1, fieldw78_1, psi78_1, evat78_1, evaw78_1, evmiut78_1, evmiuw78_1, relE78_1, conv78_1, niter78_1, mallniterc78_1, J178_1, maxgrad78_1, weight78_1] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1e-2*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1e-2*50*sech(50*(x-0.9))^2, fieldw78, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78/2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt78_1, fieldw78_1, psi78_1, evat78_1, evaw78_1, evmiut78_1, evmiuw78_1, relE78_1, conv78_1, niter78_1, mallniterc78_1, J178_1, maxgrad78_1, weight78_1] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1e-2*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1e-2*50*sech(50*(x-0.9))^2, fieldw78, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78, 1e3, 0.2, 7, 7, 1e-4, 1e3);
J190b1
J90b1
J178_1
figure
plot(t, fieldtg)
hold on
plot(t, fieldt78)
plot(t, fieldt78_1)
figure
plot(w90(1:1001), evawg(1:1001))
hold on
plot(w90(1:1001), evaw78(1:1001))
plot(w90(1:1001), evaw78_1(1:1001))
norm(psi78_1(:,end))
figure
plot(w90(1:51), fieldwg(1:51))
hold on
plot(w90(1:51), fieldw78(1:51))
plot(w90(1:51), fieldw78_1(1:51))
plot(w90(1:51), fieldw90b1(1:51))
figure
plot(t, fieldtg)
hold on
plot(t, fieldt78_1)
plot(t, fieldt90b1)
[fieldt78_2, fieldw78_2, psi78_2, evat78_2, evaw78_2, evmiut78_2, evmiuw78_2, relE78_2, conv78_2, niter78_2, mallniterc78_2, J178_2, maxgrad78_2, weight78_2] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1e-2*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1e-2*50*sech(50*(x-0.9))^2, fieldw78_1, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78_1/2, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt78_2, fieldw78_2, psi78_2, evat78_2, evaw78_2, evmiut78_2, evmiuw78_2, relE78_2, conv78_2, niter78_2, mallniterc78_2, J178_2, maxgrad78_2, weight78_2] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 1e-2*(tanh(50*(x-0.9)) - tanh(5)), @(x) 1e-2*50*sech(50*(x-0.9))^2, fieldw78_1, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
J178_2
J178_1
norm(psi78_2(:,end))
[~, ~, ~, ~, ~, ~, ~, ~, J78_2, ~, J278_2, ~, Jpnorm78_2] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 1e-2*(tanh(50*(x-0.9)) - tanh(5)), fieldw78_2, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
J278_2
J78_2
[~, ~, ~, ~, ~, ~, ~, ~, J78_2, ~, J278_2, ~, Jpnorm78_2] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 1e-2*(tanh(50*(x-0.9)) - tanh(5)), fieldw78_2, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
J78_2
J278_2
Jpnorm78_2
norm(psi78_2(:,end))
[fieldt78_2, fieldw78_2, psi78_2, evat78_2, evaw78_2, evmiut78_2, evmiuw78_2, relE78_2, conv78_2, niter78_2, mallniterc78_2, J178_2, maxgrad78_2, weight78_2] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 1e-1*(tanh(50*(x-0.9)) - tanh(5)), @(x) 1e-1*50*sech(50*(x-0.9))^2, fieldw78_1, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:51), fieldw(1:51))
plot(0:pi/T:50*pi/T, fieldw(1:51))
whos
RKlincoef
(-8-RKlincoef(2))/RKlincoef(1)
10^ans
sg8Ho = 10^((-8-3.7185e+01)/-8.7738e+00)
RK8Ho = 10^(-8-RKlincoef(2))/RKlincoef(1)
RK8Ho = 10^((-8-RKlincoef(2))/RKlincoef(1))
RK8Ho/sg8Ho
RK9Ho = 10^((-9-RKlincoef(2))/RKlincoef(1))
sg9Ho = 10^((-9-3.7185e+01)/-8.7738e+00)
RK9Ho/sg9Ho
figure
plot(0:pi/T:50*pi/T, fieldw(1:51))
sum(J1fun)
sum(J2fun)
penalnormf(normpsiT)
figure
J178_2
norm(psi78_2(:,end))
plot(w90(1:51), fieldw78_2(1:51))
figure
plot(t, fieldt78_2)
[fieldt78i, fieldw78i, psi78i, evat78i, evaw78i, evmiut78i, evmiuw78i, relE78i, conv78i, niter78i, mallniterc78i, J178i, maxgrad78i, weight78i] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) 20*(exp(-(w-0.06).^2/(2*0.005^2))).*sin((w-0.06)*pi/0.015), @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(psi78i(:,end))
norm(psig(:,end))
[fieldt78i1, fieldw78i1, psi78i1, evat78i1, evaw78i1, evmiut78i1, evmiuw78i1, relE78i1, conv78i1, niter78i1, mallniterc78i1, J178i1, maxgrad78i1, weight78i1] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw78i, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78i, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
norm(psi78i2(:,end))
norm(psi78i1(:,end))
plot(t, fieldtg)
hold on
plot(t, fieldt78i1)
figure
plot(w90(1:51), fieldwg(1:51))
hold on
plot(w90(1:51), fieldw78i1(1:51))
figure
plot(w90(1:1001), evawg(1:1001))
hold on
plot(w90(1:1001), evaw78i(1:1001))
figure
plot(w90(1:1001), log10(evawg(1:1001)^2))
plot(w90(1:1001), log10(evawg(1:1001).^2))
plot(w90(1:1001), log10(evaw78i1(1:1001).^2))
plot(w90(1:1001), log10(evawg(1:1001).^2))
hold on
plot(w90(1:1001), log10(evaw78i1(1:1001).^2))
plot(w90(1:1001), evaw78i1(1:1001))
save field78i fieldt78i1 fieldw78i1 evat78i1 evaw78i1 evmiut78i1 evmiuw78i1 relE78i1 conv78i1 niter78i1 mallniterc78i1 J178i1 maxgrad78i1 weight78i1 conv78i
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 5e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 1e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
-1.1153e-03/-2.2306e-05
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
figure
plot(t, fieldtg1)
plot(w90(1:1001), log10(evawg1(1:1001).^2))
plot(w90(1:51), fieldwg1(1:51))
norm(psig1(:,end))
plot(t, fieldtg1)
plot(t, dctI((tanh(300*(w90 - 0.02)) - tanh(200*(w90 - 0.13))).*sin((w90-0.06)*pi/0.015))/dctfactor)
hold on
plot(t, dctI((tanh(300*(w90 - 0.02)) - tanh(200*(w90 - 0.13))).*sin((w90-0.06)*pi/0.015))/dctfactor90)
plot(t, dctI((tanh(300*(w90 - 0.02)) - tanh(200*(w90 - 0.13))).*sin((w90-0.06)*pi/0.03))/dctfactor90)
plot(t, dctI((tanh(300*(w90 - 0.02)) - tanh(200*(w90 - 0.13))).*sin((w90-0.06)*pi/0.02))/dctfactor90)
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.03), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
norm(psig1(:,end))
figure
plot(t, fieldtg1)
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.03), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)])
sum([0.5*fieldw(1)*coswT(1), fieldw(2:Nt).*coswT(2:Nt), 0.5*fieldw(Nt + 1)*coswT(Nt + 1)])
sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)]) > tol
sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)])
tol
sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)]) >tol
dctfactor
(sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)])) >tol
(sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)]) >tol)
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.03), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
figure
plot(t, fieldtg1)
figure
plot(w90(1:51), fieldwg1(1:51))
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
[fieldt78e, fieldw78e, psi78e, evat78e, evaw78e, evmiut78e, evmiuw78e, relE78e, conv78e, niter78e, mallniterc78e, J178e, maxgrad78e, weight78e] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 2e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
J1g1 - 4*J2g1
J1g1 - 5*J2g1
J1g1 + 5*J2g1
J1g1 + 4*J2g1
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 8e3*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
[fieldtg1, fieldwg1, psig1, evatg1, evawg1, evmiutg1, evmiuwg1, mnitercg1, Jg1, J1g1, J2g1, Jorthg1, Jpnormg1] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-5);
Jg1, J1g1, J2g1, Jorthg1, Jpnormg1
[fieldt78e, fieldw78e, psi78e, evat78e, evaw78e, evmiut78e, evmiuw78e, relE78e, conv78e, niter78e, mallniterc78e, J178e, maxgrad78e, weight78e] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) (tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))).*sin((w-0.06)*pi/0.015), @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
norm(psi78e(:,end))
plot(t, fieldt78e)
hold on
plot(t, fieldtg1)
figure
plot(w90(1:51), fieldwg1(1:51))
hold on
plot(w90(1:51), fieldw78e(1:51))
mallniterc78e
figure
plot(w90(1:1001), evaw78e(1:1001))
hold on
plot(w90(1:1001), evawg1(1:1001))
[fieldt78e1, fieldw78e1, psi78e1, evat78e1, evaw78e1, evmiut78e1, evmiuw78e1, relE78e1, conv78e1, niter78e1, mallniterc78e1, J178e1, maxgrad78e1, weight78e1] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw78e, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78e/2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt78e1, fieldw78e1, psi78e1, evat78e1, evaw78e1, evmiut78e1, evmiuw78e1, relE78e1, conv78e1, niter78e1, mallniterc78e1, J178e1, maxgrad78e1, weight78e1] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw78e, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78e, 1e3, 0.2, 7, 7, 1e-4, 1e3);
norm(psi78e1(:,end))
figure
plot(w90(1:51), fieldw78e1(1:51))
0.9/0.1257
0.9/0.1225
relE
relE78e1
figure
plot(t, fieldtg1)
hold on
plot(t, fieldt78e1)
conv78e1(end)
[fieldt78e2, fieldw78e2, psi78e2, evat78e2, evaw78e2, evmiut78e2, evmiuw78e2, relE78e2, conv78e2, niter78e2, mallniterc78e2, J178e2, maxgrad78e2, weight78e2] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw78e1, @(w) 5e2*(tanh(300*(w - 0.02)) - tanh(200*(w - 0.13))), @(w) 10*exp(-(w-0.78).^2/(2*0.01^2)), 0, weight78e1, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
norm(psi78e2(:,end))
plot(t, fieldt78e2)
hold on
plot(t, fieldtg1)
figure