load abs_comparison2
% The Flouqet state as an initial guess:
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabs(fi0, Vf, 1, xdomain, @(x) x, @(x) 1, @(x) 0*x, fieldw0, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 2.5e-7, 2e3, 0.2, 7, 7, 1e-3);
% An arbitrary initial guess:
[fieldtb, fieldwb, psib, evmiutb, evmiuwb, relEb, convb, niterb, mallnitercb, J1b, maxgradb, weightb] = OCfalxabs(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-2*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 2e3, 0.2, 7, 7, 1e-3);