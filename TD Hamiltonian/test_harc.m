% The program was written to check the efficiency of Hillel Tal-Ezer's
% propagator, for a forced harmonic oscillator.
% You may play with the following parameters:
T = 10; Nt = 200; Nt_ts = 9; Ncheb = 9; tol=1e-5;
% Constructing the grid:
L = 8*sqrt(pi);
dx = 2*L/128;
x = (-L:dx:(L - dx)).';
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
tic
[U mniter matvecs] = TDHcheb_tsnf(@(x) x.^2/2, @(x,t) x*cos(t), [-1 188], fi0, [-L L], T, Nt, Nt_ts, Ncheb, tol);
toc
% The mean number of iterations, for a time step (should be 1, for ideal
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