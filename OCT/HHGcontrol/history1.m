load coulomb_optV240
Vf = @(x)1-1./sqrt(x.^2+1)
[~, ~, ~, ~, P240] = gsV(Vf, [-240 240], Nx240);
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
fieldwg = 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015);
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(fieldwg)/dctfactor;
figure
plot(t, fieldtg)
plot(w, fieldwg)
plot(w(1:101), fieldwg(1:101))
[~, ~, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
figure
plot(w(1:1001), evaw(1:1001))
plot(w(1:1001), evawg(1:1001))
E(1:10)-E0
whos
E240(1:10)-E0240
npsig=sqnorm(psig);
figure
plot
plot(t, npsig)
viewVPmiux(psig, Vabs240, xabs240, fieldtg, x240, 0.01)
[~, ~, psige, evatge, evawge, evmiutge, evmiuwge, mnitercge, Jge, J1ge, J2ge, Jorthge, Jpnormge] = guessresults_pnaE0b(sqrt(0.95)*P240(:,1) + sqrt(0.05)*P240(:,2), Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jge, J1ge, J2ge, Jorthge, Jpnormge]
npsige=sqnorm(psige);
figure
plot(t, npsige)
viewVPmiux(psige, Vabs240, xabs240, fieldtge, x240, 0.01)
viewVPmiux(psige, Vabs240, xabs240, fieldtg, x240, 0.01)
figure
plot(w(1:1001), evawge(1:1001))
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldte, fieldwe, psie, evate, evawe, evmiute, evmiuwe, relEe, conve, nitere, mallniterce, J1e, maxgrade, weighte] = OCfpnorm_evaE0b(sqrt(0.95)*P240(:, 1) + sqrt(0.05)*P240(:, 2), Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldte, fieldwe, psie, evate, evawe, evmiute, evmiuwe, relEe, conve, nitere, mallniterce, J1e, maxgrade, weighte] = OCfpnorm_evaE0b(sqrt(0.95)*P240(:, 1) + sqrt(0.05)*P240(:, 2), Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldte, fieldwe, psie, evate, evawe, evmiute, evmiuwe, relEe, conve, nitere, mallniterce, J1e, maxgrade, weighte] = OCfpnorm_evaE0b(sqrt(0.95)*P240(:, 1) + sqrt(0.05)*P240(:, 2), Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 5e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldte, fieldwe, psie, evate, evawe, evmiute, evmiuwe, relEe, conve, nitere, mallniterce, J1e, maxgrade, weighte] = OCfpnorm_evaE0b(sqrt(0.95)*P240(:, 1) + sqrt(0.05)*P240(:, 2), Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.77, 0.79), 0, 2e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
save excited0578strongf fieldte fieldwe evate evawe evmiute evmiuwe relEe conve nitere mallniterce J1e maxgrade weighte fieldt fieldw evat evaw evmiut evmiuw relE conv niter mallniterc J1 maxgrad weight
save excited0578strongf fieldte fieldwe evate evawe evmiute evmiuwe relEe conve nitere mallniterce J1e maxgrade weighte fieldt fieldw evat evaw evmiut evmiuw relE conv1 niter1 mallniterc J1 maxgrad weight