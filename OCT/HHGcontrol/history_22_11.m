boardw = softrectfun(w, 0.03, 0.14, 200);
plot(w(1:101), boardw(1:101))
options4 = optionsOCqn(1e-4, 1e3);
options4.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*softrectfun(w, 0.03, 0.14, 200));
[fieldtg4, fieldwg4, psig4, evatg4, evawg4, evmiutg4, evmiuwg4, mnitercg4, Jg4, J1g4, J2g4, Jorthg4, Jpnormg4] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
fieldtg4 = dctI(softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015))/dctfactor;
plot(t, fieldtg4)
[fieldtg4, fieldwg4, psig4, evatg4, evawg4, evmiutg4, evmiuwg4, mnitercg4, Jg4, J1g4, J2g4, Jorthg4, Jpnormg4] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg4, J1g4, J2g4, Jorthg4, Jpnormg4]
[fieldt4, fieldw4, psi4, evat4, evaw4, evmiut4, evmiuw4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, alpha4, invHess4] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) softrectfun(w, 0.03, 0.14, 200).*sin((w-0.085)*pi/0.015), @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options4, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:101), fieldw4(1:101))
plot(t, fieldt4)
[fieldt4a, fieldw4a, psi4a, evat4a, evaw4a, evmiut4a, evmiuw4a, relE4a, conv4a, niter4a, mallniterc4a, J14a, maxgrad4a, alpha4a, invHess4a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw4, @(w) 5e5*softrectfun(w, 0.03, 0.14, 200), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options4, 1e3, 0.2, 7, 7, 1e-4, 1e3);