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
[U1, mniter1, matvecs1, errors1, Ferror1] = TDHxp_tests(K, V, Vt, [], [-1 188], gs, x, [0 T], Nt, 5, 13, 1e-12);
plot(x, DCnewton)
plot(x, DCnewton./1i*Hpsi(K, Vts, Cnewton_der(:, 1)))
plot(x, DCnewton./(1i*Hpsi(K, Vts, Cnewton_der(:, 1))))
plot(x, (1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton)
plot(x, imag((1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton))
plot(x, (1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton)
plot(x, imag((1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton))
plot(x, (1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton)
plot(x, imag((1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton))
plot(x, (1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton)
plot(x, imag((1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton))
norm(1i*Hpsi(K, Vts, Cnewton_der(:, 1))./norm(DCnewton)
norm(1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./norm(DCnewton)
plot(x, (1i*Hpsi(K, Vts, Cnewton_der(:, 1)))./DCnewton)
figure
plot(x, DCnewton)
hold on
plot(x, (1i*Hpsi(K, Vts, Cnewton_der(:, 1))))
plot(x, norm(1i*Hpsi(K, Vts, Cnewton_der(:, 1)))/norm(DCnewton)*DCnewton)
norm(1i*Hpsi(K, Vts, Cnewton_der(:, 1)))/norm(DCnewton)
plot(x, -norm(1i*Hpsi(K, Vts, Cnewton_der(:, 1)))/norm(DCnewton)*DCnewton)
[U1, mniter1, matvecs1, errors1, Ferror1] = TDHxp_tests(K, V, Vt, [], [-1 188], gs, x, [0 T], Nt, 5, 13, 1e-12);
figure
t=0:0.05:10;
plot(t(2:end), errors1(:,10)./errors1(:,9))
plot(t(2:end), errors1(:,10)./errors1(:,8))
figure
plot(t(2:end), errors1(:,10) - errors1(:,8))
hold on
plot(t(2:end), errors1(:,9))
figure
plot(t(2:end), (errors1(:,10) - errors1(:,8))./errors1(:,9))
plot(t(2:end), errors1(:,10)./errors1(:,7))
plot(t(2:end), errors1(:,10) - errors1(:,7))
figure
plot(t(2:end), (errors1(:,10) - errors1(:,7))./errors1(:,9))
figure
plot(t, cos(t))
mx1 = evmiu(U1, x);
figure
plot(t, mx1)
Vt2 = @(u, x, t) x*cos(2*t) + conj(u).*u;
[U2, mniter2, matvecs2, errors2, Ferror2] = TDHxp_tests(K, V, Vt2, [], [-1 188], gs, x, [0 T], Nt, 5, 13, 1e-12);
mx2 = evmiu(U2, x);
figure
plot(t, mx2)
figure
plot(t(2:end), errors1(:,9))
plot(t(2:end), errors2(:,9))
plot(t(2:end), errors2(:,7))
plot(t(2:end), errors2(:,10))
plot(t(2:end), errors2(:,10)./errors2(:, 7))
figure
plot(t(2:end), (errors2(:,10)-errors2(:, 7))./errors2(:,9))
Vt3 = @(u, x, t) x*cos(0.5*t) + conj(u).*u;
[U3, mniter3, matvecs3, errors3, Ferror3] = TDHxp_tests(K, V, Vt3, [], [-1 188], gs, x, [0 T], Nt, 5, 13, 1e-12);
figure
plot(t(2:end), (errors3(:,10)-errors3(:, 7))./errors3(:,9))
figure
plot(t(2:end), errors3(:,10)./errors3(:, 7))
figure
plot(t(2:end), errors3(:,9))
Vt4 = @(u, x, t) 1.5*x*cos(t) + conj(u).*u;
[U4, mniter4, matvecs4, errors4, Ferror4] = TDHxp_tests(K, V, Vt4, [], [-1 188], gs, x, [0 T], Nt, 5, 13, 1e-12);
viewP(U4, x, 0.05)
Vt4 = @(u, x, t) 2*x*cos(t) + conj(u).*u;
[U4, mniter4, matvecs4, errors4, Ferror4] = TDHxp_tests(K, V, Vt4, [], [-1 188], gs, x, [0 T], Nt, 5, 13, 1e-12);
viewP(U4, x, 0.05)
figure
plot(t(2:end), errors4(:,9))
plot(t(2:end), errors4(:,10)./errors4(:, 7))
figure
plot(t(2:end), (errors4(:,10)-errors4(:, 7))./errors4(:,9))
save texp_error5