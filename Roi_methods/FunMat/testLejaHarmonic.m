% You may play with the following parameters:
T = 10; Nt = 100; Npcycle = 20; maxNp = 1500; tol = 1e-10;
% Constructing the grid:
Nx = 128;
L = 8*sqrt(pi);
dx = 2*L/Nx;
x = (-L:dx:(L - dx)).';
% The diagonal potential and kinetic energy operators:
V = x.^2/2;
dp = pi/L;
p = (0:dp:(2*pi/dx - dp)).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
Hoperation = @(psi) Hpsi(K, V, psi);
f = @(x) exp(-1i*x);
[U, Npf, maxEestimated] = fMtLeja(Hoperation, [0 187], f, exp(1i*8*x).*fi0, T, Nt, Npcycle, maxNp, tol);
% Computation of the maximal error - the deviation from the analytical
% result of the expectation value:
dt = T/Nt;
t=0:dt:T;
% Computing the expectation value of x in all the time points:
mx = evmiu(U, x);
figure
plot(t, mx)
xlabel('$t$ [a.u]', 'interpreter', 'latex')
ylabel('$\left<\hat\mathbf{X}\right>$ [a.u]', 'interpreter', 'latex')
error_analytic = mx - 8*sin(t);
[maxEanalytic, imaxEa] = max(abs(error_analytic))