load toda
% J1 = 9.5919e-001
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OClimf(fi0, P(:, 5), Vf, [-10 85], xdomain, @(w) sech(20*(w-1).^4), @(w) 100*sech(20*(w-1).^4), 0.1, 100, 0.05, 5, 9, 1e-3);