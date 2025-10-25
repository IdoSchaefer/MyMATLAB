% The program demonstrates the application of TDHgeneral, for a forced
% harmonic oscillator, with an inhomogeneous term.
% You may play with the following parameters:
T = 10; Nt = 400; Nt_ts = 7; Ncheb = 7; tol=1e-5;
% Constructing the grid:
L = 16*sqrt(pi);
Nx = 128;
dx = L/Nx;
x = (-L/2:dx:(L/2 - dx)).';
% Constructing the kinetic energy matrix diagonal in the p domain:
p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
Kp = p.^2/2;
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
% The time dependent disturbance operation function:
Vt = @(v, u, t, x, K) x*cos(t).*v;
ih = @(t, x, K) 0.1*cos(x*t);
tic
[U mniter matvecs] = TDHgeneral(@harmonicH0op, Vt, ih, [-1 188], fi0, [0 T], Nt, Nt_ts, Ncheb, tol, x, Kp);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which make the most of the
% numerical effort:
matvecs
% Checking if everything is OK - compairing to the ode45 result:
RKfun = @(t, u, x, K) -1i*(harmonicH0op(u, x, Kp) + Vt(u, u, t, x, Kp)) + ih(t, x, Kp);
options = odeset('RelTol', 1e-10, 'absTol', 1e-13);
tic
[time, URK] = ode45(RKfun, [0 T/2 T], fi0, options, x, Kp);
toc
URK = URK(end, :).';
error = norm(U(:, end) - URK(:, end))/norm(URK(:, end))