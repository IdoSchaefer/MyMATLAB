options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[~, ~, ~, ~, P240, ~] = gsV(Vabs240, xdomain240, Nx240);
psi0sup = sqrt(0.9)*fi0240 + sqrt(0.1)*P240(:,2);
[fieldt9, fieldw9, psi9, evat9, evaw9, evmiut9, evmiuw9, relE9, conv9, niter9, mallniterc9, J19, maxgrad9, alpha9, invHess9] = OCfx_qn(psi0sup, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
options9a = optionsOCqn(1e-4, 1e3);
options9a.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt9a, fieldw9a, psi9a, evat9a, evaw9a, evmiut9a, evmiuw9a, relE9a, conv9a, niter9a, mallniterc9a, J19a, maxgrad9a, alpha9a, invHess9a] = OCfx_qn(psi0sup, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e3*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options9a, 1e3, 0.2, 7, 7, 1e-4, 1e3);
save superposition_13
%%% Ionization:
%%% sqnorm(psi9(:,end)): 9.4308e-01
%%% sqnorm(psi9a(:,end)): 9.4601e-01
