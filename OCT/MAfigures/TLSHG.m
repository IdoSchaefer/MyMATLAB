% A TLS harmonic generation example.
miu = [0 1;
       1 0];
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfMLSrlx([1;0], [1 4], [-1 6], miu,...
    @(w) sech(20*(w-1).^4), @(w) 20*sech(20*(w-1).^4), @(w) exp(-10*(w-3).^2), 0.5, 100, 0.1, 5, 5, 1e-3);