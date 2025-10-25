options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(x) 0, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
psi(:,end).*conj(psi(:, end))
sum(psi(:,end).*conj(psi(:, end)))
figure
w
w = 0:pi/1e3:pi/0.2;
plot(w(1:101), fieldw(1:101))
t=0:0.2:1e3;
figure
plot(t, fieldt)
J1
J
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e9*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0, fieldw, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
options = optionsOCqn(1e-8, 1e3);
options = optionsOCqn(1e-4, 1e3);
options1 = optionsOCqn(1e-8, 1e3);
options1.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options1, 1e3, 0.2, 7, 7, 1e-8, 1e3
options1.invHess0 = invHess;
figure
plot(w(1:1001), evaw(1:1001))
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options1, 1e3, 0.2, 7, 7, 1e-8, 1e3);
relE1
options1.invHess0 = [];
options1 = optionsOCqn(1e-4, 1e3);
options1.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options1, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Jg, J1g, J2g, Jorthg, Jpnormg]
sum(psig(:,end).*conj(psig(:, end)))
figure
plot(t, fieldtg)
max(abs(fieldtg))
max(abs(fieldt))
sum(fieldtg)
sum(fieldtg.^2)
sum(fieldt.^2)
fieldtgn = fieldtg*sum(fieldt.^2)/sum(fieldtg.^2);
figure
plot(t, fieldtgn)
[fieldtgn, fieldwgn, psign, evatgn, evawgn, evmiutgn, evmiuwgn, mnitercgn, Jgn, J1gn, J2gn, Jorthgn, Jpnormgn] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldtgn, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
mnitercgn
sum(psign(:,end).*conj(psign(:, end)))
sum(psig(:,end).*conj(psig(:, end)))
[fieldtgn, fieldwgn, psign, evatgn, evawgn, evmiutgn, evmiuwgn, mnitercgn, Jgn, J1gn, J2gn, Jorthgn, Jpnormgn] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwgn, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[fieldtgn, fieldwgn, psign, evatgn, evawgn, evmiutgn, evmiuwgn, mnitercgn, Jgn, J1gn, J2gn, Jorthgn, Jpnormgn] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015)*sum(fieldt.^2)/sum(fieldtg.^2), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
mnitercgn
[Jgn, J1gn, J2gn, Jorthgn, Jpnormgn]
sum(psign(:,end).*conj(psign(:, end)))
plot(t, fieldtg)
hold on
plot(t, fieldt)
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3
options1 = optionsOCqn(1e-4, 1e3);
options1.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options1, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
options
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
sum(psi1(:,end).*conj(psi1(:, end)))
figure
plot(t, fieldt)
hold on
plot(t, fieldt1)
J1
J11
figure
plot(w(1:101), fieldw(1:101))
hold on
plot(w(1:101), fieldw1(1:101))
[fieldt1a, fieldw1a, psi1a, evat1a, evaw1a, evmiut1a, evmiuw1a, relE1a, conv1a, niter1a, mallniterc1a, J11a, maxgrad1a, alpha1a, invHess1a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt1b, fieldw1b, psi1b, evat1b, evaw1b, evmiut1b, evmiuw1b, relE1b, conv1b, niter1b, mallniterc1b, J11b, maxgrad1b, alpha1b, invHess1b] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw1a, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
relE1b
options
options1 = optionsOCqn(1e-5, 1e3);
options1.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
options1.invHess0 =
options1.invHess0 = invHess1b;
[fieldt1c, fieldw1c, psi1c, evat1c, evaw1c, evmiut1c, evmiuw1c, relE1c, conv1c, niter1c, mallniterc1c, J11c, maxgrad1c, alpha1c, invHess1c] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw1b, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options1, 1e3, 0.2, 7, 7, 1e-5, 1e3);
mallniterc1b
mallniterc1c
figure
plot(t, fieldt1c)
J11c
J11b
J11a
J11
J1
sum(psi1c(:,end).*conj(psi1c(:, end)))
figure
plot(w(1:1001), evaw(1:1001))
hold on
plot(w(1:1001), evaw1(1:1001))
plot(w(1:1001), evaw1c(1:1001))
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, J21c, ~, Jpnorm1c] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015)*sum(fieldt.^2)/sum(fieldtg.^2), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[J1c, J11c, J21c, Jpnorm1c]
[conv1c(end), J11c, J21c, Jpnorm1c]
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, J21c, ~, Jpnorm1c] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015)*sum(fieldt.^2)/sum(fieldtg.^2), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, J21c, ~, Jpnorm1c] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1c, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[conv1c(end), J11c, J21c, Jpnorm1c]
options1d = optionsOCqn(1e-4, 1e3);
options1c.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)));
clear options1c
options1d.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt1d, fieldw1d, psi1d, evat1d, evaw1d, evmiut1d, evmiuw1d, relE1d, conv1d, niter1d, mallniterc1d, J11d, maxgrad1d, alpha1d, invHess1d] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw1c, @(w) 4e4*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options1d, 1e3, 0.2, 7, 7, 1e-4, 1e3);
whos
figure
plot(t, fieldt1c)
iw1c = instwcos(fieldt1c, 1e3);
figure
plot(t, iw1c)
boardw = softrectfun(w, 0.03, 0.13, 50);
figure
plot(w(1:101), boardw(1:101))
boardw = softrectfun(w, 0.03, 0.13, 100);
plot(w(1:101), boardw(1:101))
boardw = softrectfun(w, 0.03, 0.14, 200);
plot(w(1:101), boardw(1:101))
[fieldt1a, fieldw1a, psi1a, evat1a, evaw1a, evmiut1a, evmiuw1a, relE1a, conv1a, niter1a, mallniterc1a, J11a, maxgrad1a, alpha1a, invHess1a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt3, fieldw3, psi3, evat3, evaw3, evmiut3, evmiuw3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, alpha3, invHess3] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt3a, fieldw3a, psi3a, evat3a, evaw3a, evmiut3a, evmiuw3a, relE3a, conv3a, niter3a, mallniterc3a, J13a, maxgrad3a, alpha3a, invHess3a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(t, fieldt3)
figure
plot(w(1:101), fieldw3(1:101))
plot(w(1:101), fieldw3a(1:101))
plot(t, fieldt3a)
sum(psi3a(:,end).*conj(psi3a(:, end)))
J13a
figure
plot(w(1:1001), evaw3a(1:1001))
[fieldt6, fieldw6, psi6, evat6, evaw6, evmiut6, evmiuw6, relE6, conv6, niter6, mallniterc6, J16, maxgrad6, alpha6, invHess6] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-1.02).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:1001), evaw6(1:1001))
sum(psi6(:,end).*conj(psi6(:, end)))
figure
plot(t, fieldt6)
whos
save fields_13_15_17
[fieldt7, fieldw7, psi7, evat7, evaw7, evmiut7, evmiuw7, relE7, conv7, niter7, mallniterc7, J17, maxgrad7, alpha7, invHess7] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.84).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt7a, fieldw7a, psi7a, evat7a, evaw7a, evmiut7a, evmiuw7a, relE7a, conv7a, niter7a, mallniterc7a, J17a, maxgrad7a, alpha7a, invHess7a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw7, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.84).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:101), fieldw7a(1:101))
figure
plot(t, fieldt7a)
figure
J17a
sum(psi7a(:,end).*conj(psi7a(:, end)))
plot(w(1:1001), evaw7a(1:1001))
iw7a = instwcos(fieldt7a, 1e3);
figure
plot(t, iw7a)
[fieldt8, fieldw8, psi8, evat8, evaw8, evmiut8, evmiuw8, relE8, conv8, niter8, mallniterc8, J18, maxgrad8, alpha8, invHess8] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.81).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt8a, fieldw8a, psi8a, evat8a, evaw8a, evmiut8a, evmiuw8a, relE8a, conv8a, niter8a, mallniterc8a, J18a, maxgrad8a, alpha8a, invHess8a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw8, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.81).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:101), fieldw8a(1:101))
figure
plot(w(1:1001), evaw8a(1:1001))
81/13
81/15
81/12
81/11
22/5
22/3
figure
plot(t, fieldt8a)
iw8a = instwcos(fieldt8a, 1e3);
figure
plot(t, iw8a)
psi0sup = sqrt(0.9)*fi0240 + sqrt(0.1)*P240(:,2);
whos
[~, ~, ~, ~, P240, ~] = gsV(Vabs240, xdomain240, Nx240);
max(abs(P240(:,1)-fi0240))
max(abs(P240(:,1)-fi0240)./abs(fi0240))
figure
plot(x, fi0240.*conj(fi0240))
plot(x240, fi0240.*conj(fi0240))
hold on
plot(x240, P240(:,1).*conj(P240(:,1)))
psi0sup = sqrt(0.9)*fi0240 + sqrt(0.1)*P240(:,2);
plot(x240, P240(:,2).*conj(P240(:,2)))
norm(psi0sup)
[fieldt9, fieldw9, psi9, evat9, evaw9, evmiut9, evmiuw9, relE9, conv9, niter9, mallniterc9, J19, maxgrad9, alpha9, invHess9] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt9, fieldw9, psi9, evat9, evaw9, evmiut9, evmiuw9, relE9, conv9, niter9, mallniterc9, J19, maxgrad9, alpha9, invHess9] = OCfx_qn(psi0sup, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
sum(psi9(:,end).*conj(psi9(:, end)))
figure
plot(w(1:1001), evaw9(1:1001))
figure
plot(t, fieldt9)
(0.073/5.33e-9)^2
(0.1/5.33e-9)^2
figure
plot(w(1:101), fieldw9(1:101))
J19
J11c
J13
J13a
options9a = optionsOCqn(1e-4, 1e3);
options9a.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt9, fieldw9, psi9, evat9, evaw9, evmiut9, evmiuw9, relE9, conv9, niter9, mallniterc9, J19, maxgrad9, alpha9, invHess9] = OCfx_qn(psi0sup, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options9a, 1e3, 0.2, 7, 7, 1e-4, 1e3);
[fieldt9a, fieldw9a, psi9a, evat9a, evaw9a, evmiut9a, evmiuw9a, relE9a, conv9a, niter9a, mallniterc9a, J19a, maxgrad9a, alpha9a, invHess9a] = OCfx_qn(psi0sup, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options9a, 1e3, 0.2, 7, 7, 1e-4, 1e3);
J19
J19a
sum(psi9a(:,end).*conj(psi9a(:, end)))
figure
plot(t, fieldt9a)
figure
plot(w(1:101), fieldw9a(1:101))
(0.03/5.33e-9)^2
figure
plot(w(1:1001), evaw9a(1:1001))
iw9a = instwcos(fieldt9a, 1e3);
figure
plot(t, iw9a)
plot(t(100:end), iw9a(100:end))
oc09a = fi0240'*psi9a;
figure
plot(t, oc09a)
plot(t, oc09a.*conj(oc09a))
oc19a = P240(:,2)'psi9a;
oc19a = P240(:,2)'*psi9a;
figure
plot(t, 0c19a)
plot(t, oc19a)
plot(t, oc19a.*conj(oc19a))
npsi9a = sqnorm(psi9a);
figure
plot(t, npsi9a)
mE9a = evH(npsi9a, @(psi) Hpsi(K240, Vabs240, psi));
mE9a = evH(psi9a, @(psi) Hpsi(K240, Vabs240, psi));
figure
plot(t, mE9a)
plot(t, mE9a./npsi9a)
whos