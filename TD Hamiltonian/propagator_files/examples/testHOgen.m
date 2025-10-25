% The program demonstrates the application of TDHgeneral, for a forced harmonic oscillator.
% You may play with the following parameters:
T = 10; Nt = 200; Nt_ts = 9; Ncheb = 9; tol=1e-5;
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
tic
[U mniter matvecs] = TDHgeneral(@harmonicH0op, Vt, [], [-1 188], fi0, [0 T], Nt, Nt_ts, Ncheb, tol, x, Kp);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which make the most of the
% numerical effort:
matvecs
% Computation of the maximal error - the deviation from the analytical
% result of the expectation value:
dt = T/Nt;
t=0:dt:T;
% Computing the expectation value of x in all the time points:
mx = evmiu(U, x);
if ~isfinite(mx(end))
    % If the result diverges, which means that we should use other
    % parameters:
    display('Error.')
end
error = mx - (-0.5*sin(t).*t);
maxer = max(abs(error))