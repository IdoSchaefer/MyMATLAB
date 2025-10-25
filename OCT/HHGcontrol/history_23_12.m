doc dct
Vf = @(x) 1 - 1./sqrt(1+x.^2);
function fundamentalw = get_fundamental(Vf, fieldv, xdomain, Nx, m)
40/64
function fundamentalw = get_fundamental(Vf, fieldv, xdomain, Nx, m)
fundamentalw = get_fundamental(Vf, 0:0.001:0.1, [-200, 200], 640);
0.01^2/3.55e-9
0.1^2/3.55e-9
0.1^2/5.33e-9
(0.01/5.33e-9)^2
(0.001/5.33e-9)^2
(0.05/5.33e-9)^2
figure
plot(0:0.001:0.1, fundamentalw)
fundamentalw = get_fundamental(Vf, 0:0.01:0.1, [-200, 200], 640);
plot(0:0.01:0.1, fundamentalw)
fundamentalw = get_fundamental(Vf, 0:0.001:0.1, [-200, 200], 640);
plot(0:0.001:0.1, fundamentalw)
plot((0:0.001:0.1).^2/5.33e-9, fundamentalw)
plot(((0:0.001:0.1)/5.33e-9).^2, fundamentalw)
plot(log10(((0:0.001:0.1)/5.33e-9).^2), fundamentalw)
fundamentalw = get_fundamental(Vf, sqrt(10^(10.5:0.1:14))*5.33e-9, [-200, 200], 640);
fundamentalw = get_fundamental(Vf, sqrt(10.^(10.5:0.1:14))*5.33e-9, [-200, 200], 640);
plot(10.5:0.1:14, fundamentalw)
sqrt(10^11.6)*5.33e-9
(0.01/5.33e-9)^2
3.363*200
3.363e-3*200
whos
fundamentalw = get_fundamental(Vabs240, sqrt(10.^(10.5:0.1:14))*5.33e-9, [-240, 240], 768);
plot(10.5:0.1:14, fundamentalw)
function fundamentalw = get_fundamental(Vf, fieldv, xdomain, Nx, m)
[fi0, E0, x, E, P, H] = gsV(Vabs240, [-240 240], 768);
fundamentalw = get_fundamental1(Vabs240, xabs240, sqrt(10.^(10.5:0.1:14))*5.33e-9, [-240, 240]);
plot(10.5:0.1:14, fundamentalw)
[fi0, E0, x, E, P, H] = gsV(Vabs240-3.3630e-03*xabs240, [-240 240], 768);
figure
plot(x, conj(fi0).*fi0)
plot(x, conj(P(:, 2)).*P(:,2))
plot(x, conj(P(:, 3)).*P(:,3))
plot(x, conj(P(:, 4)).*P(:,4))
plot(x, conj(P(:, 5)).*P(:,5))
plot(x, conj(P(:, 6)).*P(:,6))
plot(x, conj(P(:, 7)).*P(:,7))
plot(x, conj(P(:, 8)).*P(:,8))
plot(x, conj(P(:, 9)).*P(:,9))
plot(x, conj(P(:, 10)).*P(:,10))
plot(x, conj(P(:, 11)).*P(:,11))
plot(x, conj(P(:, 12)).*P(:,12))
sqrt(10^11.5)*5.33e-9
[fi0, E0, x, E, P, H] = gsV(Vabs240-2.9973e-03*xabs240, [-240 240], 768);
plot(x, conj(fi0).*fi0)
plot(x, conj(P(:, 2)).*P(:,2))
sqrt(10^11.7)*5.33e-9
[fi0, E0, x, E, P, H] = gsV(Vabs240-3.7734e-03*xabs240, [-240 240], 768);
plot(x, conj(fi0).*fi0)
plot(x, conj(P(:, 2)).*P(:,2))
plot(x, conj(P(:, 3)).*P(:,3))
plot(x, conj(P(:, 4)).*P(:,4))
plot(x, conj(P(:, 5)).*P(:,5))
plot(x, conj(P(:, 6)).*P(:,6))
plot(x, conj(P(:, 7)).*P(:,7))
plot(x, conj(P(:, 8)).*P(:,8))
3.363e-3*40
whos
fundamentalw = get_fundamental1(Vabs, xabs, sqrt(10.^(10.5:0.1:14))*5.33e-9, [-80 80]);
plot(10.5:0.1:14, fundamentalw)
sqrt(10^12.2)*5.33e-9
ans*40
sqrt(10^12.3)*5.33e-9
ans*40
E0
sqrt(10^13.1)*5.33e-9
sqrt(10^13)*5.33e-9
ans*40
1.8912e-02*40
figure
plot(x, conj(P(:, 8)).*P(:,8))
plot(x, conj(P(:, 1)).*P(:,1))
plot(x, conj(P(:, 2)).*P(:,2))
[fi0, E0, x, E, P, H] = gsV(Vabs-1.8912e-2*xabs, [-80 80], 768);
[fi0, E0, x, E, P, H] = gsV(Vabs-1.8912e-2*xabs, [-80 80], 256);
plot(x, conj(P(:, 1)).*P(:,1))
plot(x, conj(P(:, 2)).*P(:,2))
plot(x, conj(P(:, 3)).*P(:,3))
plot(x, conj(P(:, 4)).*P(:,4))
plot(x, conj(P(:, 5)).*P(:,5))
plot(x, conj(P(:, 6)).*P(:,6))
plot(x, conj(P(:, 7)).*P(:,7))
plot(x, conj(P(:, 8)).*P(:,8))
plot(x, conj(P(:, 9)).*P(:,9))
plot(x, conj(P(:, 10)).*P(:,10))
plot(x, conj(P(:, 11)).*P(:,11))
plot(x, conj(P(:, 12)).*P(:,12))
plot(x, conj(P(:, 13)).*P(:,13))
plot(x, conj(P(:, 14)).*P(:,14))
plot(x, conj(P(:, 15)).*P(:,15))
plot(x, conj(P(:, 16)).*P(:,16))
plot(x, conj(P(:, 17)).*P(:,17))
plot(x, conj(P(:, 18)).*P(:,18))
E(1:18)
0.703-0.329
whos
clear all
options = optionsOCqn(1e-3, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-4, 1e3
figure
dctfactor = 1e3/(sqrt(5e3*pi))
w = 0:pi/1e3:pi/0.2;
t = 0:0.2:1e3;
fieldgw = 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015);
fieldg1w = 5*exp(-(w-0.06).^2/(2*0.01^2)).*cos((w-0.06)*pi/0.015);
plot(w(1:101), fieldgw(1:101))
fieldgt = dctI(fieldgw)/dctfactor;
fieldg1t = dctI(fieldg1w)/dctfactor;
figure
plot(t, fieldgt)
hold on
plot(t, fieldg1t)
plot(t, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)))/dctfactor;
plot(t, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)))/dctfactor);
(0.05/5.33e-9)^2
(0.048/5.33e-9)^2
figure
plot(w, exp((w-0.78).^2/(2*0.005^2))
plot(w, exp((w-0.78).^2/(2*0.005^2)))
plot(w, exp(-(w-0.78).^2/(2*0.005^2)))
dw = pi/1e3;
0.78/dw
plot(w(200:300), exp(-(w(200:300)-0.78).^2/(2*0.005^2)))
plot(w(200:300), exp(-(w(200:300)-0.78).^4/(2*0.005^2)))
plot(w(200:300), exp(-(w(200:300)-0.78).^4/(2*0.001^2)))
plot(w(200:300), exp(-(w(200:300)-0.78).^4/(2*0.0005^2)))
plot(w(200:300), exp(-(w(200:300)-0.78).^2/(2*0.005^2)))
sum(exp(-(w(200:300)-0.78).^2/(2*0.005^2)))*dw
sum(rectanglefun(w, 0.77, 0.79))*dw
sum(exp(-(w(200:300)-0.78).^2/(2*0.003^2)))*dw
sum(exp(-(w(200:300)-0.78).^2/(2*0.007^2)))*dw
sum(exp(-(w(200:300)-0.78).^2/(2*0.008^2)))*dw
plot(w(200:300), exp(-(w(200:300)-0.78).^2/(2*0.008^2)))
0.01/13
plot(w(200:300), exp(-(w(200:300)-0.78).^2/(2*0.01^2)))
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)));
pi/0.015
plot(t, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(w*200))/dctfactor);
plot(t, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*200))/dctfactor);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
options
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) 0, @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 4e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e9*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e9*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e9*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
min_interval
b-a
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 9e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e9*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
boardw = softrectfun(w, 0.03, 0.14, 200);
plot(w(1:101), boardw(1:101))
options4 = optionsOCqn(1e-4, 1e3);
options4.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*softrectfun(w, 0.03, 0.14, 200);
options4.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*softrectfun(w, 0.03, 0.14, 200));
[fieldtg4, fieldwg4, psig4, evatg4, evawg4, evmiutg4, evmiuwg4, mnitercg4, Jg4, J1g4, J2g4, Jorthg4, Jpnormg4] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) @(w) 5e5*softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[fieldtg4, fieldwg4, psig4, evatg4, evawg4, evmiutg4, evmiuwg4, mnitercg4, Jg4, J1g4, J2g4, Jorthg4, Jpnormg4] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) @(w) 5e5*softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[fieldtg4, fieldwg4, psig4, evatg4, evawg4, evmiutg4, evmiuwg4, mnitercg4, Jg4, J1g4, J2g4, Jorthg4, Jpnormg4] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
fieldtg4 = dctI(softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015))/dctfactor;
plot(t, fieldtg4)
[fieldtg4, fieldwg4, psig4, evatg4, evawg4, evmiutg4, evmiuwg4, mnitercg4, Jg4, J1g4, J2g4, Jorthg4, Jpnormg4] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg4, J1g4, J2g4, Jorthg4, Jpnormg4]
[fieldt3, fieldw3, psi3, evat3, evaw3, evmiut3, evmiuw3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, alpha3, invHess3] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options4, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt4, fieldw4, psi4, evat4, evaw4, evmiut4, evmiuw4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, alpha4, invHess4] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options4, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:101), fieldw4(1:101))
plot(t, fieldt4)
[fieldt4a, fieldw4a, psi4a, evat4a, evaw4a, evmiut4a, evmiuw4a, relE4a, conv4a, niter4a, mallniterc4a, J14a, maxgrad4a, alpha4a, invHess4a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw4, @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options4, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt5, fieldw5, psi5, evat5, evaw5, evmiut5, evmiuw5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*(exp(-(w-0.06).^2/(2*0.01^2)) + 0.1*(exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.01^2)) + exp(-(w-0.18).^2/(2*0.01^2)))), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options5, 1e3, 0.2, 7, 7, 1e-4, 1e3);
sum(psi5(:,end).*conj(psi5(:,end)))
sum(psi4a(:,end).*conj(psi4a(:,end)))
J14a
J15
figure
plot(t, fieldt5)
figure
plot(w(1:101), fieldw5(1:101))
figure
plot(t, fieldt5)
plot(t, fieldt4a)
figure
plot(w(1:101), fieldw4a(1:101))
78/5
78/4
pi/1e3*sqrt(2/pi)
pi/1e3*sqrt(2/pi)/2
sqrt(2/5e3)/dctfactor
(0.01/5.33e-9)^2
(0.01*1.2533e-03/5.33e-9)^2
(0.56*1.2533e-03/5.33e-9)^2
whos
hold on
plot(w(1:101), broadw(1:101))
plot(w(1:101), boardw(1:101))
iw4a = instwcos(fieldt4a, 1e3);
figure
plot(t, iw4a)
78/4
iw5 = instwcos(fieldt5, 1e3);
figure
plot(t, iw5)
fieldw5fft = 1/(sqrt(2*pi))*0.2*fft(fieldt5(1:end-1));
fieldw4afft = 1/(sqrt(2*pi))*0.2*fft(fieldt4a(1:end-1));
figure
plot(-pi/0.2:2*pi/1e3:(pi/0.2 - 2*pi/1e3), fieldw4afft)
plot(0:2*pi/1e3:50*2*pi/1e3, fieldw4afft(1:51)
plot(0:2*pi/1e3:50*2*pi/1e3, fieldw4afft(1:51))
figure
plot(w(1:101), fieldw4a(1:101))
plot(w(1:101), fieldw4a(1:101).^2)
hold on
plot(0:2*pi/1e3:50*2*pi/1e3, fieldw4afft(1:51).*conj(fieldw4afft(1:51)))
plot(0:2*pi/1e3:50*2*pi/1e3, fieldw4afft(1:51).*conj(fieldw4afft(1:51))*6)
max(fieldw4a(1:101).^2)
max(fieldw4a(1:101).^2)/max(fieldw4afft(1:51).*conj(fieldw4afft(1:51)))
plot(0:2*pi/1e3:50*2*pi/1e3, fieldw4afft(1:51).*conj(fieldw4afft(1:51))*4.8761e+00)
npsi4a = sqnorm(psi4a);
figure
plot(t, npsi4a)
whos
gsoc4a = fi0240'*psi4a;
figure
plot(t, gsoc4a)
gsoc4a = (fi0240'*psi4a).*conj(fi0240'*psi4a);
plot(t, gsoc4a)
colors1 = exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2));
colors1t = dctI(colors1)/dctfactor;
figure
plot(w(1:101), colors1(1:101))
options5a = optionsOCqn(1e-4, 1e3);
options5a.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2))));
[fieldt5a, fieldw5a, psi5a, evat5a, evaw5a, evmiut5a, evmiuw5a, relE5a, conv5a, niter5a, mallniterc5a, J15a, maxgrad5a, alpha5a, invHess5a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2))), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options5a, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt5a, fieldw5a, psi5a, evat5a, evaw5a, evmiut5a, evmiuw5a, relE5a, conv5a, niter5a, mallniterc5a, J15a, maxgrad5a, alpha5a, invHess5a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.005^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2))), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options5a, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldtg4a, fieldwg4a, psig4a, evatg4a, evawg4a, evmiutg4a, evmiuwg4a, mnitercg4a, Jg4a, J1g4a, J2g4a, Jorthg4a, Jpnormg4a] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 10*exp(-(w-0.06).^2/(2*0.005^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg4a, J1g4a, J2g4a, Jorthg4a, Jpnormg4a]
figure
plot(t, fieldtg4a)
[fieldt5a, fieldw5a, psi5a, evat5a, evaw5a, evmiut5a, evmiuw5a, relE5a, conv5a, niter5a, mallniterc5a, J15a, maxgrad5a, alpha5a, invHess5a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 10*exp(-(w-0.06).^2/(2*0.005^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2))), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options5a, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
J15a
sum(psi5a(:,end).*conj(psi5a(:,end)))
plot(w(1:101), fieldw5a(1:101))
figure
plot(t, fieldt5a)
iw5a = instwcos(fieldt5a, 1e3);
figure
plot(t, iw5a)
options5a = optionsOCqn(1e-4, 1e3);
options5a.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2))));
options5b = optionsOCqn(1e-4, 1e3);
options5b.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2)) + exp(-(w-0.18).^2/(2*0.005^2))));
[fieldt5b, fieldw5b, psi5b, evat5b, evaw5b, evmiut5b, evmiuw5b, relE5b, conv5b, niter5b, mallniterc5b, J15b, maxgrad5b, alpha5b, invHess5b] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 10*exp(-(w-0.06).^2/(2*0.005^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2)) + exp(-(w-0.18).^2/(2*0.005^2))), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options5b, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt5c, fieldw5c, psi5c, evat5c, evaw5c, evmiut5c, evmiuw5c, relE5c, conv5c, niter5c, mallniterc5c, J15c, maxgrad5c, alpha5c, invHess5c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw5b, @(w) 5e5*(exp(-(w-0.06).^2/(2*0.005^2)) + exp(-(w-0.03).^2/(2*0.005^2)) + exp(-(w-0.12).^2/(2*0.005^2)) + exp(-(w-0.18).^2/(2*0.005^2))), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options5b, 1e3, 0.2, 7, 7, 1e-4, 1e3);
conv4b(end)
J14b
figure
plot(w(1:101), fieldw4b(1:101))
figure
plot(t, fieldt4b)
figure
plot(w(1:1001), evaw4b(1:1001))
figure
plot(w(1:101), fieldw5b(1:101))
figure
plot(t, fieldt4b)
plot(t, fieldt5b)
plot(t, fieldt5c)
plot(w(1:101), fieldw5c(1:101))
sum(psi5c(:,end).*conj(psi5c(:,end)))
ocgs5c = fi0240'*psi5c;
figure
plot(t, ocgs5c)
ocgs5c = (fi0240'*psi5c).*conj((fi0240'*psi5c));
plot(t, ocgs5c)
P240
[~, ~, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
Vf = @(x) 1 - 1./sqrt(1+x.^2);
[~, ~, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
ocfund5c = (P240(:,2)'*psi5c).*conj((P240(:,2)'*psi5c));
figure
plot(t, ocfund5c)
figure
plot(t, ocgs5c+ocfund5c)
W = psi2wigner(fieldt5c(1:10:end));
W = psi2wigner(fieldt5c(1:10:end).');
clear W
W5c = psi2wigner(fieldt5c(1:10:end).');
contourf(W5c)
contourf(W5c(500:600, :))
contourf(W5c(500:700, :))
contourf(W5c(500:600, :))
contourf(0:1:1e3, 0:2*pi/1e3:100*2*pi/1e3,W5c(500:600, :))
contourf(0:2*pi/1e3:100*2*pi/1e3, 0:1:1e3, W5c(500:600, :))
contourf(0:2*pi/1e3:100*2*pi/1e3, 0:1:1e3, W5c(501:600, :))
W5c = psi2wigner(fieldt5c(1:10:end-1).');
contourf(W5c)
contourf((0:2*pi/1e3:100*pi/1e3).', 0:1:1e3, W5c(501:550, :))
contourf((0:2*pi/1e3:100*pi/1e3).', 0:1:999, W5c(501:550, :))
contourf((0:2*pi/1e3:98*pi/1e3).', 0:1:999, W5c(501:550, :))
contourf((0:2*pi/1e3:98*pi/1e3), 0:1:999, W5c(501:550, :))
contourf(0:1:999, (0:2*pi/1e3:98*pi/1e3), W5c(501:550, :))
size(W5c)
size(fieldt5c)
figure
contourf(W5c)
figure
sum(fieldt7a)*0.2
sum(fieldt7a.^2)*0.2
plot(t, 0.1*exp(-(t-500).^2./(2*100^2))
plot(t, 0.1*exp(-(t-500).^2./(2*100^2)))
plot(t, 0.1*exp(-(t-500).^2./(2*100^2)).*sin(0.06*t))
plot(t, 0.1*exp(-(t-500).^2./(2*100^2)).*sin(0.06*(t-500)))
plot(t, 0.1*exp(-(t-500).^2./(2*200^2)).*cos(0.06*(t-500)))
plot(t, 0.1*exp(-(t-500).^2./(2*150^2)).*cos(0.06*(t-500)))
sum(0.1*exp(-(t-500).^2./(2*150^2)).*cos(0.06*(t-500)))*0.2
sum((0.1*exp(-(t-500).^2./(2*150^2)).*cos(0.06*(t-500))).^2)*0.2
sum((0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500))).^2)*0.2
plot(t, 0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500)))
plot(t, 0.1*exp(-(t-500).^2./(2*150^2)).*cos(0.06*(t-500)))
xlabel('$t$ (a.u.)', 'interpreter', 'latex')
ylabel('$\epsilon(t)$ (a.u.)', 'interpreter', 'latex')
whos
figrue
figure
size(conv1d)
size(conv1c)
size(conv1b)
size(conv1a)
conv1a
conv1b(end)
conv1(end)
conv1b(1)
conv(1c)
conv1c
plot(w(1:1001), evaw1b(1:1001))
sum((0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500))).^2)
sum(fieldt1b.^2)
figure
plot(t, fieldt1b)
max(abs(fieldt1b))
trans_lim = 0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500));
figure
plot(t, trans_lim)
max(abs(trans_lim))
figure
trans_limw = dctI(trans_lim)*dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi));
trans_limw = dctI(trans_lim)*dctfactor;
figure
figure
plot(w(1:101), trans_lim(1:101))
plot(w(1:101), trans_limw(1:101))
dw = pi/1e3;
trans_limw_con = fieldw20b(trans_limw, 5e5*exp(-(w-0.06).^2/(2*0.01^2)), dw);
trans_limw_con = fieldw20b(trans_limw, 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw);
trans_limw_con = fieldw20b(trans_limw.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw);
figure
plot(w(1:101), trans_limw_con(1:101))
hold on
plot(w(1:101), trans_limw(1:101))
figure
trans_lim_con = dctI(trans_limw_con)/dctfactor;
plot(t, trans_lim_con)
hold on
plot(t, trans_lim)
sum(fieldt3a.^2)
sum(fieldt1b.^2)
sum(fieldt6.^2)
sum(fieldt7a.^2)
sum((0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500))).^2)
max(abs(trans_lim))
max(abs(fieldt1b))
max(abs(fieldt3a))
max(abs(fieldt6))
max(abs(fieldt7a))
[fieldttl, fieldwtl, psitl, evattl, evawtl, evmiuttl, evmiuwtl, mniterctl, Jtl, J1tl, J2tl, ~, Jpnormtl] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), trans_limw_con*max(abs(fieldt1b))/max(abs(trans_lim_con)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[fieldttl, fieldwtl, psitl, evattl, evawtl, evmiuttl, evmiuwtl, mniterctl, Jtl, J1tl, J2tl, ~, Jpnormtl] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), trans_limw_con.'*max(abs(fieldt1b))/max(abs(trans_lim_con)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
max(abs(fieldt1b))/max(abs(trans_lim_con))
figure
plot(t, fieldt1b)
hold on
plot(t, fieldttl)
Jtl
J1tl
figure
plot(w(1:101), fieldwtl(1:101))
Jpnormtl
npsitl = sqnorm(psitl);
figure
plot(t, npsitl)
figure
plot(w(1:1001), evawtl(1:1001))
hold on
plot(w(1:1001), evaw1b(1:1001))
[fieldttl, fieldwtl, psitl, evattl, evawtl, evmiuttl, evmiuwtl, mniterctl, Jtl, J1tl, J2tl, ~, Jpnormtl] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), trans_limw_con.', @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
J1tl
Jpnormtl
figure
plot(w(1:1001), evawtl(1:1001))
npsitl = sqnorm(psitl);
figure
plot(t, npsitl)
figure
plot(t, fieldttl)
hold on
plot(t, fieldttl)
plot(w(1:1001), evawtl(1:1001))
plot(w(1:1001), evaw3a(1:1001))
plot(w(1:1001), evaw6(1:1001))
plot(w(1:1001), evaw7a(1:1001))
figure
plot(w(1:1001), evaw7a(1:1001))
dw
0.06/dw
19*dw
5.97*13
5.97*15
5.97*17
trans_lim1w = 10*exp(-(w-0.06).^2/(2*dw^2)).*sin((w-0.06)*500);
trans_lim1 = dctI(trans_lim1w)/dctfactor;
figure
plot(t, trans_lim1)
cw = 0.1*cos(0.06*(t-500));
plot(t, cw)
cww = dctI(cw);
cww = dctI(cw)*dctfactor;
plot(w(1:1001), cww(1:1001))
cw = 0.1*sin(0.06*(t-500));
cww = dctI(cw)*dctfactor;
plot(w(1:1001), cww(1:1001))
cwwfil = cww*exp(-(w-0.06).^2/(2*0.01^2));
cwwfil = cww.*exp(-(w-0.06).^2/(2*0.01^2));
plot(w(1:1001), cwwfil(1:1001))
plot(w(1:101), cwwfil(1:101))
cwfil = dctI(cwwfil)/dctfactor;
plot(t, cwfil)
max(abs(cwfil))
cww_con.' = fieldw20b(cwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw);
cww_con.' = fieldw20b(cwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw).';
cww_con = fieldw20b(cwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw).';
figure
plot(w(1:101), cw_con)
plot(w(1:101), cww_con(1:101))
cww_con = fieldw20b(cwwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw).';
plot(w(1:101), cww_con(1:101))
cw_con = dctI(cww_con)/dctfactor;
plot(t, cw_con)
2*pi/0.06
[fieldtcw, fieldwcw, psicw, evatcw, evawcw, evmiutcw, evmiuwcw, mniterccw, Jcw, J1cw, J2cw, ~, Jpnormcw] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
Jcw
J1cw
figure
sum(cw_con.^2)
plot(w(1:1001), evawcw(1:1001))
npsicw = sqnorm(psitl);
figure
plot(t, npsicw)
hold on
plot(w(1:1001), evaw1b(1:1001))
plot(w(1:1001), evawcw(1:1001))
hold on
plot(w(1:1001), evaw1b(1:1001))
plot(w(1:1001), evaw3a(1:1001))
plot(w(1:1001), evaw6(1:1001))
sum((0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500))).^2)
sum(fieldt3a.^2)
max(abs(fieldt1b))/max(abs(cw_con))
[fieldtcw1b, fieldwcw1b, psicw1b, evatcw1b, evawcw1b, evmiutcw1b, evmiuwcw1b, mniterccw1b, Jcw1b, J1cw1b, J2cw1b, ~, Jpnormcw1b] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*max(abs(fieldt1b))/max(abs(cw_con)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
J1cw1b
figure
plot(w(1:1001), evawcw1b(1:1001))
hold on
plot(w(1:1001), evaw1b(1:1001))
sum(fieldtcw1b.^2)
max(abs(fieldtcw1b))
[fieldtcw6, fieldwcw6, psicw6, evatcw6, evawcw6, evmiutcw6, evmiuwcw6, mniterccw6, Jcw6, J1cw6, J2cw6, ~, Jpnormcw6] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*max(abs(fieldt6))/max(abs(cw_con)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-1.02).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
max(abs(fieldtcw6))
figure
plot(t, fieldt6)
plot(w(1:1001), evawcw6(1:1001))
J1cw6
hold on
plot(w(1:1001), evaw6(1:1001))
J16
[fieldtcw3a, fieldwcw3a, psicw3a, evatcw3a, evawcw3a, evmiutcw3a, evmiuwcw3a, mniterccw3a, Jcw3a, J1cw3a, J2cw3a, ~, Jpnormcw3a] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*max(abs(fieldt3a))/max(abs(cw_con)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
J1cw3a
J13a
figure
plot(w(1:1001), evawcw3a(1:1001))
hold on
plot(w(1:1001), evaw3a(1:1001))
sum(fieldtcw3a.^2)
sum(fieldt3a.^2)
sum(fieldtcw3a.^2)/sum(fieldt3a.^2)
sum(fieldt3a.^2)/sum(fieldtcw3a.^2)
sqrt(sum(fieldt3a.^2)/sum(fieldtcw3a.^2))
[fieldtcwe3a, fieldwcwe3a, psicwe3a, evatcwe3a, evawcwe3a, evmiutcwe3a, evmiuwcwe3a, mniterccwe3a, Jcwe3a, J1cwe3a, J2cwe3a, ~, Jpnormcwe3a] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*sqrt(sum(fieldt3a.^2)/sum(cw_con.^2)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
plot(w(1:1001), evawcwe3a(1:1001))
Jpnormcwe3a
sum(fieldt3a.^2)
sum(fieldtcwe3a.^2)
figure
plot(t, fieldtcwe6)
plot(t, fieldtcwe3a)
0.05^2
(0.05/5.33e-9)^2
J13a
J1cwe3a
max(abs(fieldt1b))
max(abs(fieldt3a))
max(abs(fieldt6))
figure
plot(t, fieldt1b)
plot(t, fieldt3a)
plot(t, fieldt6)
plot(t, fieldt7a)
[fieldtcwe1b, fieldwcwe1b, psicwe1b, evatcwe1b, evawcwe1b, evmiutcwe1b, evmiuwcwe1b, mniterccwe1b, Jcwe1b, J1cwe1b, J2cwe1b, ~, Jpnormcwe1b] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_con*sqrt(sum(fieldt1b.^2)/sum(cw_con.^2)), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
npsicw1b = sqnorm(psicw1b);
figure
plot(t, npsicw1b)
plot(t, npsicwe1b)
npsicwe1b = sqnorm(psicw1b);
plot(t, npsicwe1b)
npsicwe1b = sqnorm(psicwe1b);
plot(t, npsicwe1b)
J1cw1b
J1cwe1b
J11b
whos
clear psicw psicw1b psicw3a psicw6 psicwe1b psicwe3a psig psign psitl psi1c psi1d
save fields_13_15_17