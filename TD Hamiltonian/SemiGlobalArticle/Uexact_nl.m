[fi0240b, E0240b, x240b, E240b, P240b, H240b] = gsV(Vabs240, xdomain240, Nx240);
[gs, niter] = gsNLHdiag(H240b, @(u,x,t) conj(u).*u, x240, 1e-12);
[Uex, mniterex, matvecsex] = TDHxpKr1(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech(-(t-500)/(170)).^2.*cos(0.06*(t-500)) + conj(u).*u, [], gs, x240, 0:1000, 3e4, 9, 13, eps);