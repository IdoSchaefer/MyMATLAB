% The program assumes the existence of the variables:
% T, Nt, Nt_ts, Ncheb, tol, gs.
T = 10; Nt = 200; Nt_ts = 9; Ncheb = 9; tol=1e-5;
% Constructing the grid:
L = 16*sqrt(pi);
Nx = 128;
dx = L/Nx;
x = (-L/2:dx:(L/2 - dx)).';
tic
[U mniter matvecs] = NLHchebn(@(x) x.^2/2, @(u, x, t) conj(u).*u + x*cos(t), [-1 188], gs, [-L/2 L/2], T, Nt, Nt_ts, Ncheb, tol);
%[U mniter matvecs] = NLHchebn(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-1 188], exp(1i*x).*gs, [-L/2 L/2], T, Nt, Nt_ts, Ncheb, tol);
toc
mniter
matvecs
dt = T/Nt;
t=0:dt:T;
%mx = evx(U, x);
mx = evmiu(U, x);
if ~isfinite(mx(end))
    display('Error.')
end
%error = (mx - sin(t));
error = mx - (-0.5*sin(t).*t);
max(abs(error))