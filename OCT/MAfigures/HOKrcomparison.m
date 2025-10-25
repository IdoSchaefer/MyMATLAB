xdomain = [-8*sqrt(pi) 8*sqrt(pi)];
dx = 16*sqrt(pi)/128;
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
tar = pi^(-1/4)*exp(-(x - 1).^2/2)*sqrt(dx);
[field, psi, relE, conv, niter, convprop, npropconv, mallniterc, J1, maxgrad, weight] = OCrlxcp(fi0, tar, @(x) x.^2/2,...
    [-10 188], xdomain, @(t) 1, 0.1, 1e-1, 10, 0.05, 5, 9, 1e-3);
[fieldK, psiK, relEK, convK, niterK, mniterK] = OCcheb(fi0, tar, @(x) x.^2/2, [-10 188], xdomain, @(t) 1, 0.1, 10, 0.05, 5, 9, 1e-3);