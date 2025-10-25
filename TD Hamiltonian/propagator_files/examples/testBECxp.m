% The program demonstrates the application of TDHxp,
% for a nonlinear Hamiltonian - BEC harmonic trap, under a forcing perturbation.
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
K = p.^2/2;
% The potential energy matrix diagonal in the x domain:
V = x.^2/2;
% We have to construct the H0 matrix, to find the ground state (it's unnecessary
% for the propagation process).
% The potential energy matrix:
Vmat = diag(V);
% The kinetic energy matrix in the x domain:
Kmat = Nx*ifft(ifft(diag(K))')';
% The Hamiltonian:
H = Kmat + Vmat;
% The ground state, found by an iterative process:
gs = gsNLHdiag(H, @(u,x,t) conj(u).*u, x, 2e-12);
% The time dependent disturbance operation function:
Vt = @(u, x, t) x*cos(t) + conj(u).*u;
tic
[U mniter matvecs] = TDHxp(K, V, Vt, [], [-1 188], gs, x, [0 T], Nt, Nt_ts, Ncheb, tol);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which make the most of the
% numerical effort:
matvecs
% Computation of the maximal error - the deviation from the expected
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