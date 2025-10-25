xdomain = [-8*sqrt(pi) 8*sqrt(pi)];
dx = 16*sqrt(pi)/128;
x = (xdomain(1):dx:xdomain(2)-dx).';
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
tar = (4/pi)^(1/4)*x.*exp(-x.^2/2)*sqrt(dx);
[field, psi, relE, conv, niter, mniter] = OCcheb(fi0, tar, @(x) x.^2/2, [-10 188], xdomain, @(t) t.^0, 0.05, 20, 0.05, 5, 9, 1e-3);