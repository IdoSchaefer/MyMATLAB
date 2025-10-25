load HCl02_2
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1785, [-0.08 0.25], xdomain, miuf4, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2), @(w) 100*rectanglefun(w, 3.6e-2, 3.8e-2), 10, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-2);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1785, [-0.08 0.25], xdomain, miuf4, P(:, 1:20), P(:, 21:32), fieldw, @(w) 2500*rectanglefun(w, 0, 1.5e-2), @(w) 100*rectanglefun(w, 3.6e-2, 3.8e-2), 10, @(n) n.^2, weight, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv2, niter2, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1785, [-0.08 0.25], xdomain, miuf4, P(:, 1:20), P(:, 21:32), fieldw, @(w) 2500*rectanglefun(w, 0, 1.5e-2), @(w) 100*rectanglefun(w, 3.6e-2, 3.8e-2), 10, @(n) n.^2, weight, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv3, niter3, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1785, [-0.08 0.25], xdomain, miuf4, P(:, 1:20), P(:, 21:32), fieldw, @(w) 2500*rectanglefun(w, 0, 1.5e-2), @(w) 100*rectanglefun(w, 3.6e-2, 3.8e-2), 10, @(n) n.^2, weight, 10000, 5, 5, 9, 1e-3);
nitertot=niter+niter1+niter2+niter3;
convtot=[conv conv1(2:end) conv2(2:end) conv3(2:end)];
% J1= 2.6731e-001
% Df0val = 1.9260e-001
% miuf4 = @(x) 0.5*Df0val*x.*(1-tanh(x-0.7)) = 9.6300e-002*x.*(1-tanh(x-0.7))