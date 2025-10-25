% A 3LS harmonic generation example.
miu = [0 1 1;
       1 0 1;
       1 1 0];
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfMLSrlx([1;0;0], [1 1.9 3], [-1 5], miu,...
    @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), 1, 100, 0.05, 5, 5, 1e-2);
%    @(w) rectanglefun(w, 0.85, 1.15), @(w) 0.25*rectanglefun(w, 0.85, 1.15), @(w) rectanglefun(w, 1.85, 2.15), 1, 100, 0.05, 5, 5, 1e-2);
% J1 = 2.1502e+001
