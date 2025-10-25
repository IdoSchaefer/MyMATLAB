% The program assumes the existence of the variables:
% T, Nt, Nt_ts, tol, gs.
dx = 16/128;
x = (-8:dx:7.875).';
% The potential energy matrix:
V = diag(x.^2/2);
dp = 2*pi/16;
p = (-pi/dx):dp:(pi/dx - dp);
% The kinetic energy matrix in the p domain:
Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
K = 128*ifft(ifft(ifftshift(Kp))')';
% The Hamiltonian:
H = K + V;
% % The harmonic oscillator ground state vector:
% fi0 = pi^(-1/4)*exp(-x.^2/2)/sqrt(8);
tic
[U mniter] = NLHdiagn(H, @(u, x, t) conj(u).*u, exp(1i*x).*gs, x, T, Nt, Nt_ts, tol);
toc
mniter
matvecs = mniter*Nt*3 + Nt + Nt_ts
dt = T/Nt;
t=0:dt:T;
mx = evx(U, x);
error = mx - sin(t);
max(abs(error))