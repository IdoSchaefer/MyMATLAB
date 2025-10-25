%function [gs, niter] = gsBEC(H0, Vnl, x, tol)
% dx = 16/128;
% x = (-8:dx:7.875).';
dx = sqrt(2*pi/128);
xdlength = 128*dx;
x = (-xdlength/2:dx:(xdlength/2 - dx)).';
% The potential energy matrix:
V = diag(x.^2/2);
% dp = 2*pi/16;
dp = 2*pi/xdlength;
p = (-pi/dx):dp:(pi/dx - dp);
% The kinetic energy matrix in the p domain:
Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
K = 128*ifft(ifft(ifftshift(Kp))')';
% The Hamiltonian:
H0 = K + V;
[gs, niter] = gsNLHdiag(H0, @(u, x, t)  conj(u).*u, x, 2e-12);