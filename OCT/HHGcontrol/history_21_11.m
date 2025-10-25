options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
w = 0:pi/1e3:pi/0.2;
t=0:0.2:1e3;
plot(w(1:101), fieldw(1:101))
plot(t, fieldt)
figure
[fieldtg2, fieldwg2, psig2, evatg2, evawg2, evmiutg2, evmiuwg2, mnitercg2, Jg2, J1g2, J2g2, Jorthg2, Jpnormg2] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, fieldw, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg2, J1g2, J2g2, Jorthg2, Jpnormg2]
options2 = optionsOCqn(1e-4, 1e3);
options2.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options2, 1e3, 0.2, 7, 7, 1e-4, 1e3);
hold on
fieldwg = 5e5*exp(-(w-0.06).^2/(2*0.01^2));
fieldtg = dctI(fieldwg)/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(fieldwg)/dctfactor;
plot(t, fieldt)
hold on
plot(t, fieldtg)
fieldwg = 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015);
fieldtg = dctI(fieldwg)/dctfactor;
plot(t, fieldtg)
[fieldtg2, fieldwg2, psig2, evatg2, evawg2, evmiutg2, evmiuwg2, mnitercg2, Jg2, J1g2, J2g2, Jorthg2, Jpnormg2] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, fieldw/2, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg2, J1g2, J2g2, Jorthg2, Jpnormg2]
options2 = optionsOCqn(1e-4, 1e3);
options2.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw/2, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options2, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt2a, fieldw2a, psi2a, evat2a, evaw2a, evmiut2a, evmiuw2a, relE2a, conv2a, niter2a, mallniterc2a, J12a, maxgrad2a, alpha2a, invHess2a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options2, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt2b, fieldw2b, psi2b, evat2b, evaw2b, evmiut2b, evmiuw2b, relE2b, conv2b, niter2b, mallniterc2b, J12b, maxgrad2b, alpha2b, invHess2b] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2a, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options2, 1e3, 0.2, 7, 7, 1e-4, 1e3);
J12b
J12b-conv2b(end)
figure
plot(w(1:101), fieldw1b(1:101))
plot(w(1:101), fieldw2b(1:101))
figure
plot(t, fieldt2b)
figure
plot(w(1:1001), evaw2b(1:101))
plot(w(1:1001), evaw2b(1:1001))
hold on
filterE = exp(-(w-0.06).^2/(2*0.01^2));
plot(w(1:101), filterE(1:101)*20.78)
[fieldt2c, fieldw2c, psi2c, evat2c, evaw2c, evmiut2c, evmiuw2c, relE2c, conv2c, niter2c, mallniterc2c, J12c, maxgrad2c, alpha2c, invHess2c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2b, @(w) 2e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
options3 = optionsOCqn(1e-4, 1e3);
options3.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt2c, fieldw2c, psi2c, evat2c, evaw2c, evmiut2c, evmiuw2c, relE2c, conv2c, niter2c, mallniterc2c, J12c, maxgrad2c, alpha2c, invHess2c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2b, @(w) 2e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
h
options3.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 2e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt2c, fieldw2c, psi2c, evat2c, evaw2c, evmiut2c, evmiuw2c, relE2c, conv2c, niter2c, mallniterc2c, J12c, maxgrad2c, alpha2c, invHess2c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2b, @(w) 2e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
options3.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 2e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt2c, fieldw2c, psi2c, evat2c, evaw2c, evmiut2c, evmiuw2c, relE2c, conv2c, niter2c, mallniterc2c, J12c, maxgrad2c, alpha2c, invHess2c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2b, @(w) 8e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[Jg2, J1g2, J2g2, Jorthg2, Jpnormg2]
J22b = J12b-conv2b(end)
conv2b(end)/J22b
5.2/4
4/5.2
4/5.24
options3.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 7.7e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt2c, fieldw2c, psi2c, evat2c, evaw2c, evmiut2c, evmiuw2c, relE2c, conv2c, niter2c, mallniterc2c, J12c, maxgrad2c, alpha2c, invHess2c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2b, @(w) 7.7e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
conv2b(end)-7.7*J22b
conv2b(end)-5.2*J22b
4/5.2
4e4/5.2
J12b/J22b
4/6.24
4/6.2
options3.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 6.45e3*exp(-(w-0.06).^2/(2*0.01^2)));
J12b-6.45*J22b
J12b-0.645*J22b
J12b-1/0.645*J22b
J12b-6.45/4*J22b
J12b-6.45*J22b
J12b-6.24*J22b
J12b-6.2*J22b
[fieldt2c, fieldw2c, psi2c, evat2c, evaw2c, evmiut2c, evmiuw2c, relE2c, conv2c, niter2c, mallniterc2c, J12c, maxgrad2c, alpha2c, invHess2c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, fieldw2b, @(w) 6.45e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
sum(psi2b(:,end).*conj(psi2b(:,end))
sum(psi2b(:,end).*conj(psi2b(:,end)))
0.01^2/5.33e-9
(0.01/5.33e-9)^2
(0.02/5.33e-9)^2
(0.05/5.33e-9)^2
(0.08/5.33e-9)^2
mE = evHt(psi2b, K240, Vabs240, fieldt2b, xabs240);
figure
plot(t, mE)
mE0 = evHt(psi2b, K240, Vabs240, zeros(5001,1), xabs240);
hold on
plot(t, mE0)
whos