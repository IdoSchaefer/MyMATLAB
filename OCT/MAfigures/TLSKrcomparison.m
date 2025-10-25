% The Ncheb is not ideal, mniterc>2.
[field, psi, relE, conv, niter, convprop, npropconv, mallniterc, J1, maxgrad, weight] = OCMLSrlxcp([1;0], [0;1], [1 2], [-1 4], miu,...
    @(t) 1, 0.1, 10, 0.1, 5, 3, 1e-3);
[fieldc, psic, relEc, convc, niterc, mniterc] = OCchebMLS([1;0], [0;1], [1 2], [-1 4], miu, @(t) 1, 0.1, 10, 0.1, 5, 3, 1e-4);