% You may play with the following parameters:
T = 100; Nt = 1000; Npcycle = 20; maxNp = 1500; tol = 1e-13;
Vf = @(x) 12.5*(exp(-0.2*x) - 1).^2;
% Constructing the grid:
Nx = 32;
[xdomain, dx, x] = nlVgrid(Vf, Nx);
xdlength = xdomain(2) - xdomain(1);
% The diagonal potential and kinetic energy operators:
V = Vf(x.');
dp = 2*pi/xdlength;
p = (0:dp:(2*pi/dx - dp)).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
% The ground state:
[fi0, E0, x, E, P, H] = gsV(Vf, xdomain, Nx);
Eanalytic = 12.5 - (25 - 0.5 - (0:31).').^2*0.2^2/2;
errorE = abs((E - Eanalytic)./Eanalytic);
figure
plot(0:9, log10(errorE(1:10)))
xlabel('$n$', 'interpreter', 'latex')
ylabel('$\log_{10}(\epsilon_{rel})$', 'interpreter', 'latex')
Hoperation = @(psi) Hpsi(K, V, psi);
f = @(x) exp(-1i*x);
[U, Npf, maxEestimated] = fMtLeja(Hoperation, [0 22], f, exp(1i*0.5*x).*fi0, T, Nt, Npcycle, maxNp, tol);
dt = T/Nt;
t=0:dt:T;
% Computing the expectation value of x in all the time points:
mx = evmiu(U, x);
figure
plot(t, mx)
xlabel('$t$ [a.u]', 'interpreter', 'latex')
ylabel('$\left<\hat\mathbf{X}\right>$ [a.u]', 'interpreter', 'latex')
% %windowfun = 0.5*(tanh(t-5)-tanh(t-95));
% dw = pi/T;
% w = 0:dw:pi/dt;
% %ftmx = fft(mx(1:end-1).*windowfun(1:end-1));
% dctmx = dctI(mx-mean(mx));
% figure
% plot(w, dctmx)
% %plot(w, abs(ftmx(1:end/2)))
% Using the Krylov Basis Diagonalization Method to extract the
% characteristic Bohr frequencies of the system from the signal:
c = mx.';
Kre = 250;
FBDMrinp
Ndisp = 15;
[d_resulto, orderd] = sort(d_result);
w_od = wreal(orderd);
Nd = length(d_result);
wdisp = w_od((Nd-Ndisp+1):Nd);
ddisp = d_resulto((Nd-Ndisp+1):Nd);
[wdisp, orderwd] = sort(wdisp);
ddisp = ddisp(orderwd);
[wdisp.', ddisp.']
figure
bar(wdisp, ddisp)
xlabel('$\omega$ [a.u.]', 'interpreter', 'latex')
ylabel('$|A|$ [a.u]', 'interpreter', 'latex')
w1 = E(2:10) - E(1:9);
w2 = E(3:10) - E(1:8);
w3 = E(4:10) - E(1:7);
w4 = E(5:10) - E(1:6);
reler1 = abs((wdisp(5:-1:3).' - w1(1:3))./w1(1:3))
reler2 = abs((wdisp(10:-1:6).' - w2(1:5))./w2(1:5))
reler3 = abs((wdisp(13:-1:11).' - w3(1:3))./w3(1:3))
reler4 = abs((wdisp(15:-1:14).' - w4(1:2))./w4(1:2))

% reler1 =
% 
%   1.2548e-006
%   2.8838e-004
%   4.8646e-003
% 
% reler2 =
% 
%   6.1639e-007
%   1.2929e-005
%   3.7693e-004
%   7.4070e-003
%   9.2135e-004
% 
% reler3 =
% 
%   5.3368e-007
%   4.8719e-006
%   6.0794e-005
% 
% reler4 =
% 
%   1.5497e-007
%   1.3658e-006
