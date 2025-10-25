load toda
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1, [-10 85], xdomain, miuexp, P(:, 1:88), P(:, 89:128), @(w) 5*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 4.7, 4.9), 1, @(n) n, 1e-2, 100, 0.05, 9, 9, 1e-2);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1, [-10 85], xdomain, miuexp, P(:, 1:88), P(:, 89:128), fieldw, @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 4.7, 4.9), 1, @(n) n, weight, 100, 0.05, 9, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv2, niter2, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1, [-10 85], xdomain, miuexp, P(:, 1:88), P(:, 89:128), fieldw, @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 4.7, 4.9), 1, @(n) n, weight, 100, 0.05, 9, 9, 1e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv3, niter3, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1, [-10 85], xdomain, miuexp, P(:, 1:88), P(:, 89:128), fieldw, @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 4.7, 4.9), 1, @(n) n, weight, 100, 0.05, 9, 9, 1e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv4, niter4, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1, [-10 85], xdomain, miuexp, P(:, 1:88), P(:, 89:128), fieldw, @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 4.7, 4.9), 1, @(n) n, weight, 100, 0.05, 9, 9, 1e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv5, niter5, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1, [-10 85], xdomain, miuexp, P(:, 1:88), P(:, 89:128), fieldw, @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 4.7, 4.9), 1, @(n) n, weight, 100, 0.05, 9, 9, 1e-4);
convtot=[conv conv1(2:end) conv2(2:end) conv3(2:end) conv4(2:end) conv5(2:end)];
nitertot=niter+niter1+niter2+niter3+niter4+niter5;
% nitertot=1224
% relE = 1.7760e-004

