% The program assumes the existence of the variables:
% T, Nt, gs.
dx = sqrt(2*pi/128);
xdlength = 128*dx;
x = (-xdlength/2:dx:(xdlength/2 - dx)).';
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p(65:128) = p(65:128) - 2*pi/dx;
K = p.^2/2;
dt = T/Nt;
%fi0 = pi^(-1/4)*exp(-x.^2/2)/sqrt(8);
tic
uf = RK4uf(@SEnl, [0 T], exp(1i*8*x).*gs, dt, K, @(u, x, t) x.^2/2 + conj(u).*u, x);
%U = RK4(@SEnl, [0 T], gs, dt, K, @(u, x, t) x.^2/2 + conj(u).*u + x*cos(t), x);
toc
matvecs = Nt*4;
mx = evx(uf, x);
if ~isfinite(mx(end))
    display('Error.')
end
%error = (mx - sin(t));
error = (mx - 8*sin(T));
%error = mx - (-0.5*sin(t).*t);
abs(error)