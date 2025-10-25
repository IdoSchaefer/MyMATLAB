[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalx(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-6, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfalx(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), fieldw, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, weight, 2e3, 0.2, 5, 9, 5e-4);
covtot = [conv conv1(2:end)];
nitertot = niter + niter1;
[fieldtb, fieldwb, psib, evmiutb, evmiuwb, relEb, convb, niterb, mallnitercb, J1b, maxgradb, weightb] = OCfalx(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 3e-2*(-0.5*(tanh(x+35)-tanh(x-35))+1), fieldw, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, weight, 2e3, 0.2, 5, 9, 1e-4);
[fieldtb, fieldwb, psib, evmiutb, evmiuwb, relEb, convb1, niterb1, mallnitercb, J1b, maxgradb, weightb] = OCfalx(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 3e-2*(-0.5*(tanh(x+35)-tanh(x-35))+1), fieldwb, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, weightb, 2e3, 0.2, 5, 9, 1e-4);
nitertotb = niterb+niterb1;
convtotb=[convb convb1(2:end)];


