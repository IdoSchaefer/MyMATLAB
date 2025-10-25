T = 10; Nt = 1e4; Ncheb = 7;
dt = T/Nt;
tmiddle = dt/2:dt:(T - dt/2);
t=0:dt:T;
% Constructing the grid:
L = 16*sqrt(pi);
Nx = 128;
dx = L/Nx;
x = (-L/2:dx:(L/2 - dx)).';
% Constructing the kinetic energy matrix diagonal in the p domain:
p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
% The potential energy matrix diagonal in the x domain:
V = x.^2/2;
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
Hop = @(v, t) Hpsi(K, V + x*cos(t), v);
U = SchrPWCcheb(Hop, fi0, tmiddle, [-1, 188], T, Nt, Ncheb);
mx = evmiu(U, x);
error = mx - (-0.5*sin(t).*t);
maxer = max(abs(error))