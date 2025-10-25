options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt6, fieldw6, psi6, evat6, evaw6, evmiut6, evmiuw6, relE6, conv6, niter6, mallniterc6, J16, maxgrad6, alpha6, invHess6] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-1.02).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
%%% sqnorm(psi3a(:, end)): 9.3664e-01

