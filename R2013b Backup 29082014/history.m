%-- Unknown date --%
old on
plot(-0.25:0.25:10.25, imag(Voptx2), 'r')
[allT2, allR2, allnorm2] = allTRcoef1(Voptx2,  [-0.25 10.25], [1 31], 301);
figure
plot(1:0.1:31, allT2, 'r')
plot(1:0.1:31, allR2)
doc spline
save
%-- 27/11/2013 13:39 --%
load
whos
clear all
defaultop = optimset('fmincon');
optionsgrad = optimset(defaultopcon, 'Display', 'off', 'TolFun', tol, 'TolX', tol, 'Algorithm', 'trust-region-reflective', 'GradObj','on');
defaultopcon = optimset('fmincon');
optionsgrad = optimset(defaultopcon, 'Display', 'off', 'TolFun', tol, 'TolX', tol, 'Algorithm', 'trust-region-reflective', 'GradObj','on');
defaultop
optionsgrad = optimset(defaultopcon, 'Display', 'off', 'Algorithm', 'trust-region-reflective', 'GradObj','on');
[Vopt, optval] = fmincon(@(V) Vabs_efficiency1(@(V) percosV_grad(V, [0 10], [1 31], 61), -rand(42,1), [], [], [], [], [], [], [], optionsgrad)
[Vopt, optval] = fmincon(@(V) percosV_grad(V, [0 10], [1 31], 61), -rand(42,1), [], [], [], [], [], [], [], optionsgrad)
doc fminunc
options = optimoptions(@fminunc,'GradObj','on')
[Vopt, optval] = fminunc(@(V) percosV_grad(V, [0 10], [1 31], 61), -rand(42,1), options)
[Vopt, optval] = fminunc(@(V) percosV_grad(V, [0 10], [1 31], 61, 1e-5, 1, 100), -rand(42,1), options)
[performance, gradient] = percosV_grad(-rand(42,1),[0 10], [1 31], 61, 1e-5, 1, 100);
performance
[performance, gradient] = percosV_grad(-rand(42,1),[0 10], [1 31], 61, 1e-5, 1, 100);
Vxreim
format short e
Vxreim
[performance, gradient] = percosV_grad(-rand(42,1),[0 10], [1 31], 61, 1e-5, 1, 100)
[performance, gradient] = percosV_grad(dctI(-rand(42,1)),[0 10], [1 31], 61, 1e-5, 1, 100)
[performance, gradient] = percosV_grad(dctI(-ones(42,1)),[0 10], [1 31], 61, 1e-5, 1, 100)
[performance, gradient] = percosV_grad(dctI(-rand(42,1)),[0 10], [1 31], 61, 1e-5, 1, 100)
Vrand = -rand(42,1);
figure
[performance, gradient] = percosV_grad([dctI(-rand(21,1)); dctI(-rand(21,1))),[0 10], [1 31], 61, 1e-5, 1, 100)
[performance, gradient] = percosV_grad([dctI(-rand(21,1)); dctI(-rand(21,1))],[0 10], [1 31], 61, 1e-5, 1, 100)
[Vopt, optval] = fminunc(@(V) percosV_grad(V, [0 10], [1 31], 61, 1e-5, 1, 100), [dctI(-rand(21,1)); dctI(-rand(21,1))], options)
Voptx = dctI(Vopt(1:21)) + 1i*dctI(Vopt(22:42));
figure
plot(0:0.5:10, real(Voptx))
plot(0:0.5:10, imag(Voptx))
plot(0:0.5:10, real(Voptx))
[parameter, gradient] = Vabs_efficiency1(Voptx, xdomain, kdomain, Nk)
[parameter, gradient] = Vabs_efficiency1(Voptx, [0 10], [1 31], 61)
1e-5*sum(Vopt.^2)
imag(Voptx)
exp(ans*100)
sum(ans)
[Vopt, optval, exitflag, output] = fminunc(@(V) percosV_grad(V, [0 5], [1 11], 11, 1e-5, 1, 100), [dctI(-rand(11,1)); dctI(-rand(11,1))], options)
Voptx = dctI(Vopt(1:11)) + 1i*dctI(Vopt(21:22));
Voptx = dctI(Vopt(1:11)) + 1i*dctI(Vopt(12:22));
figure
plot(0:10, real(Voptx))
plot(0:10, imag(Voptx))
[Vopt, optval, exitflag, output] = fminunc(@(V) percosV_grad(V, [0 5], [1 11], 11, 1e-7, 1, 100), [dctI(-rand(11,1)); dctI(-rand(11,1))], options)
Voptx = dctI(Vopt(1:11)) + 1i*dctI(Vopt(12:22));
figure
plot(0:10, real(Voptx))
plot(0:10, imag(Voptx))
[Vopt, optval, exitflag, output] = fminunc(@(V) percosV_grad(V, [0 5], [1 11], 11, 1e-10, 1, 100), [dctI(-rand(11,1)); dctI(-rand(11,1))], options)
Voptx = dctI(Vopt(1:11)) + 1i*dctI(Vopt(12:22));
plot(0:10, real(Voptx))
plot(0:10, imag(Voptx))
options = optimoptions(@fminunc,'GradObj','on', 'OutputFcn', @getconv)
[Vopt, optval, exitflag, output] = fminunc(@(V) percosV_grad(V, [0 5], [1 11], 11, 1e-10, 1, 100), [dctI(-rand(11,1)); dctI(-rand(11,1))], options)
options = optimoptions(@fminunc,'GradObj','on', 'OutputFcn', 'PlotFcns',{@optimplotfval,@optimplotfirstorderopt});
options = optimoptions(@fminunc,'GradObj','on', 'OutputFcn', 'PlotFcns',@optimplotfval);
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt});
[Vopt, optval, exitflag, output] = fminunc(@(V) percosV_grad(V, [0 5], [1 11], 11, 1e-10, 1, 100), [dctI(-rand(11,1)); dctI(-rand(11,1))], options)
[performance, gradient] = percosV_grad([dctI(-rand(21,1)); dctI(-rand(21,1))],[0 10], [1 31], 61, 1e-5, 1, 100)
grad_norms
pospenal*t*exp(Vximag*t)
[Vximag pospenal*t*exp(Vximag*t)]
gradient
[Vopt, optval, exitflag, output] = fminunc(@(V) percosV_grad(V, [0 5], [1 11], 11, 0, 0, 100), [dctI(-rand(11,1)); dctI(-rand(11,1))], options)
[performance, gradient] = percosV_grad([dctI(-rand(21,1)); dctI(-rand(21,1))],[0 10], [1 31], 61, 0, 0, 100)
gradient
(Nx + 2):(2*Nx - 1)
grad_norms
Vxreal + 1i*Vximag
[parameter, gradient] = Vabs_efficiency1(-rand(21, 1) - 1i*rand(21,1), [0 10], [1 31], 61)
[parameter, gradient] = Vabs_efficiency1(-ones(21, 1) - 1i*ones(21,1), [0 10], [1 31], 61)
[Vopt, optval, exitflag, output] = fminunc(@(V) Vabs_efficiency1(V(1:end/2) + 1i*V(end/2+1:end), [0 5], [1 11], 11), -rand(22,1), options)
allgradnorm(:, 1:7)
gradient
[parameter, gradient] = Vabs_efficiency1(-rand(21, 1) - 1i*rand(21,1), [0 10], [1 31], 61)
[performance, gradient] = percosV_grad([dctI(-rand(11,1)); dctI(-rand(11,1))],[0 10], [1 11], 11, 0, 0, 100)
grad_perf_Vx-grad_norms
[grad_norms [2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)]]
[grad_norms, [2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)]]
[grad_norms(1:Nx), [2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)]]
dctI([2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)])
[dctI(ans), [2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)]]
[grad_norms(Nx+1:2*Nx), [2*grad_perf_Vx(Nx + 1); grad_perf_Vx((Nx + 2):(2*Nx - 1)); 2*grad_perf_Vx(2*Nx)]]
%-- 30/11/2013 18:36 --%
prefdir
dir C:\Users\Ido\AppData\Roaming\MathWorks\MATLAB\R2013a
%-- 17/12/2013 09:55 --%
prefdir
0.0009
%-- 17/12/2013 20:09 --%
path
%-- 17/12/2013 20:23 --%
%-- 18/12/2013 11:56 --%
load
whos
meshgrid(W1)
mesh(W1)
mesh(W2)
mesh(real(W2))
sum(sum(W2))
sum(sum(real(W2)))
W2 = psi2wigner(fi0);
norm(psi)
norm(psi_int)
W
mesh(real(W2))
mesh(real(W))
W2 = psi2wigner(fi0);
figure
mesh(real(W))
size(psi)
mesh(real(W(1:(2*N - 2), :)))
figure
[maxW, imax] = max(max(W))
[maxW, imax] = max(W)
fftW=fft(W)/(2*N);
W=fftshift(W);
fftW=fftshift(W);
fftW=fft(W)/(2*N);
fftW=fftshift(fftW);
mesh(real(fftW))
W
max(imag(W2))
max(max(imag(W2)))
max(max(abs(imag(W2))))
W2 = psi2wigner(fi0);
figure
mesh(real(ro))
figure
mesh(real(W))
W
W = [W(end,:); W(1:end-1, :)]
fft(W)
whos
sum(sum(ans))
W2(end, :)
%-- 19/12/2013 10:25 --%
load
W
v
ro
qi = ones(2*N, 1)*(1:2*N);
si = (1:2*N).'*ones(1, 2*N);
%     xi = qi + si - N;
%     xtagi = qi - si + N;
xi = qi + si - 1 - N;
xtagi = qi - si + 1 + N;
is_in_ro = xi>=1 & xi<=2*N & xtagi>=1 & xtagi<=2*N;
ro_v_indices = (xi - 1)*(2*N) + xtagi;
W0 = zeros(2*N);
W(is_in_ro) = ro(ro_v_indices(is_in_ro));
W0(is_in_ro) = ro(ro_v_indices(is_in_ro));
W
W0
load
W
W0
fftW0 = fft(W0)
sum(fftW0)
whos
W5 = psi2wigner(fi0);
figure
mesh(real(W5))
W5 = psi2wigner(fi0);
W(1, :)
W5 = psi2wigner(fi0);
mesh(real(W5))
max(max(abs(imag(W5))))
sum(sum(W5))
sum(W5,2)
max(abs(ans))
sumW5=sum(W5,2);
[maxsW5, imax] = max(abs(sum(W5)))
[maxsW5, imax] = max(abs(sumW5))
sum(sumW5)
W5 = psi2wigner(fi0);
figure
mesh(W)
W5 = psi2wigner(fi0);
mesh(W5)
mesh(real(W5))
W5 = psi2wigner(fi0);
mesh(real(W5))
sum(sum(W5))
W0 = psi2wigner(v);
W0
sum(sum(W))
sum(v.^2)
sum(sum(W5))
fi0p = 1/sqrt(256)*fft(fftshift(fi0));
figure
dx=16/256
x=-8:dx:(8-dx);
plot(x, fi0.*conj(fi0))
size(x)
size(fi0)
dx=16/128
x=-8:dx:(8-dx);
fi0p = 1/sqrt(128)*fft(fftshift(fi0));
plot(x, fi0.*conj(fi0))
hold on
plot(x, fi0p.*conj(fi0p), 'r')
fi0p = fftshift(fi0p);
plot(x, fi0p.*conj(fi0p), 'r')
norm(fi0)
norm(fi0p)
p=-pi/dx:pi/16:(pi/dx - pi/16);
plot(p, fi0p.*conj(fi0p), 'r')
size(p)
size(x)
p=-pi/dx:2*pi/16:(pi/dx - pi/16);
plot(p, fi0p.*conj(fi0p), 'r')
plot(p, fi0p.*conj(fi0p)*dp, 'r')
dp=pi/16;
plot(x, fi0.*conj(fi0)*dx)
plot(p, fi0p.*conj(fi0p)*dp, 'r')
plot(x, fi0.*conj(fi0)/dx)
plot(p, fi0p.*conj(fi0p)/dp, 'r')
dp=2*pi/16;
plot(p, fi0p.*conj(fi0p)/dp, 'r')
max(abs(fi0.*conj(fi0)*dx - fi0p.*conj(fi0p)/dp))
figure
dx
dp
fi0pW = sum(W5, 2);
max(abs(fi0p-fi0pW))
max(abs(fi0p-fi0pW(1:2:end)))
figure
plot(-2*pi/(dx):dp/2:(2*pi/dx - dp/2))
plot(-2*pi/(dx):dp/2:(2*pi/dx - dp/2), fi0pW)
plot(-2*pi/(dx):dp:(2*pi/dx - dp/2), fi0pW)
hold on
plot(p, fi0p.*conj(fi0p), 'r')
plot(p, fi0p.*conj(fi0p)/2, 'r')
plot(-2*pi/(dx):dp:(2*pi/dx - dp), fi0pW)
hold on
plot(p, fi0p.*conj(fi0p)/2, 'r')
plot(p*2, fi0p.*conj(fi0p)/2, 'g')
clf
hold on
plot(p, fi0p.*conj(fi0p)/2, 'r')
plot(-pi/(dx):dp:(pi/dx - dp), fi0pW)
plot(-pi/dx:dp/2:(pi/dx - dp/2), fi0pW)
fi0W=sum(W5);
figure
plot(x, fi0.*conj(fi0)/2)
hold on
plot(x, W5, 'r')
plot(-8:dx/2:(8-dx/2), W5, 'r')
clf
plot(x, fi0.*conj(fi0)/2)
hold on
plot(-8:dx/2:(8-dx/2), fi0W, 'r')
xW = -8:dx/2:(8-dx/2);
pW = -pi/dx:dp/2:(pi/dx-dp/2);
mesh(xW, pW, W5)
mesh(xW, pW, real(W5))
axis equal
W5 = psi2wigner(fi0);
L = 8*sqrt(pi);
dx = 2*L/128;
x = (-L:dx:(L - dx)).';
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
W5 = psi2wigner(fi0);
mesh(xW, pW, real(W5))
xW = -L:dx/2:(L-dx/2);
dp=pi/L;
pW = -pi/dx:dp/2:(pi/dx-dp/2);
mesh(xW, pW, real(W5))
figure
mesh(xW, pW, real(W5))
mesh(xW, pW, W5)
W1 = psi2wigner(exp(1i*x)*fi0);
W1 = psi2wigner(exp(1i*x).*fi0);
figure
mesh(xW, pW, W1)
W1 = psi2wigner(exp(1i*x).*fi0);
mesh(xW, pW, W1)
sum(sum(W1))
W5 = psi2wigner(fi0);
mesh(xW, pW, W5)
sum(sum(W5))
save
%-- 23/12/2013 11:58 --%
load
whos
Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
K = 128*ifft(ifft(ifftshift(Kp))')';
% The Hamiltonian:
H = K + V;
V = x.^2/2;
H = K + V;
V = diag(x.^2/2);
H = K + V;
[E, P] = diag(H);
[E, P] = eig(H);
[E, orderE]=sort(E)
[P, D] = eig(H);
E = diag(D);
%    [E, orderE] = sort(E);
[Ereal, orderE] = sort(real(E));
E = E(orderE);
P = P(:, orderE);
E0 = E(1);
E
size(x)
p
E/E0
E0*2
p=-pi/dx:pi/L:(pi/dx - pi/L);
V = diag(x.^2/2);
dp = 2*pi/16;
p = (-pi/dx):dp:(pi/dx - dp);
% The kinetic energy matrix in the p domain:
Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
K = 128*ifft(ifft(ifftshift(Kp))')';
% The Hamiltonian:
H = K + V;
size(V)
size(K)
dx
dp=pi/L;
p = (-pi/dx):dp:(pi/dx - dp);
% The kinetic energy matrix in the p domain:
Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
K = 128*ifft(ifft(ifftshift(Kp))')';
% The Hamiltonian:
H = K + V;
[P, D] = eig(H);
E = diag(D);
%    [E, orderE] = sort(E);
[Ereal, orderE] = sort(real(E));
E = E(orderE);
P = P(:, orderE);
E0 = E(1);
E0
real(E(1:10))
W6 = psi2wigner(P(:,1));
mesh(xW, pW, W6)
W6 = psi2wigner(P(:,2));
mesh(xW, pW, W6)
W6 = psi2wigner(P(:,3));
mesh(xW, pW, W6)
pW = -pi/dx:dp/2:(pi/dx-dp/2);
viewPp(exp(1i*x).*fi0, dx, 2*L)
%-- 25/12/2013 15:17 --%
whos
xdomain
dx
clear all
Vf
Vf1 = @(x)1-1./sqrt(x.^2+1)
whos
[U mniter matvecs] = TDHcheb_tsnf(Vf1(x), @(x,t) -0.1*x*sin(t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 5, 9, 1e-5);
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 5, 9, 1e-5);
mniter
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5);
mniter
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(t), [-2, 15.5], fi0, xdomain, 1e3, 1e4, 5, 5, 1e-5);
mniter
E(2)-E(1)
[U mniter matvecs] = TDHcheb_tsnf(Vf, @(x,t) -0.1*x*sin(t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 5, 9, 1e-5);
mniter
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(0.06*t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 5, 9, 1e-5);
mniter
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(0.06*t), [-2, 15.5], fi0, xdomain, 1e3, 1e4, 5, 9, 1e-5);
mniter
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(0.06*t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5);
mniter
viewP(U, dx, 0.01)
size(fi0)
dx=160/256
viewP(U, dx, 0.01)
viewP(U, x, 0.01)
Wrange = getWrange(U(:,1:20:end))
20*0.2
viewWigner(U, xdomain, 0.01, Wrange)
viewWigner(U(:,1:20:end), xdomain, 0.01, Wrange)
W0 = psi2wigner(fi0);
figure
surf(W, 'LineStyle', none)
surf(W0, 'LineStyle', none)
surf(W0, 'LineStyle', '')
surf(W0, 'LineStyle')
surf(W0)
surf(W0, 'LineStyle', '.')
surf(W0, 'LineStyle', 'none')
viewWigner(U(:,1:20:end), xdomain, 0.01, Wrange)
Wrange
Wrange = getWrange(U(:,1:20:end))
max(max(W))
min(min(W))
maxW = max(W);
[maxW, imaxc] = max(W);
imaxc
maxW
viewWigner(U(:,1:20:end), xdomain, 0.01, Wrange)
maxW = max(W);
maxW
Wrange = getWrange(U(:,1:20:end))
max(max(W))
viewWigner(U(:,1:20:end), xdomain, 0.01, Wrange)
min_x = xdomain(1);
max_x = xdomain(2);
xdlength = max_x - min_x;
dx = xdlength/Nx;
x = (min_x:dx:(max_x - dx)).';
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
[U mniter matvecs] = TDHxpKr(K, Vf, Vtfun, ihfun, ui, x, tdomain, Nts, Nt_ts, Nkr, tol, varargin
[U mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(0.06*t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5
[U mniter matvecs] = TDHxpKr(K, Vf, @(x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U mniter matvecs] = TDHxpKr(K, Vf, @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U mniter matvecs] = TDHxpKr(K, Vf(x), @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[Ua mniter matvecs] = TDHxpKr(K, Vf(x), @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
Wrangea = getWrange(Ua(:,1:20:end))
viewWigner(Ua(:,1:20:end), xdomain, 0.01, Wrangea)
viewVP(Ua(:,1:20:end), Vf, 0.1*sin(0:4:1e3), x)
true(2,1)
M=ones(10);
M(:,3:7) = 0;
is_none0 = M~=0
tic,sum(sum(M)),toc
tic,sum(sum(M));,toc
tic,sum(sum(M(is_none0))),toc
tic,sum(sum(M(is_none0)));,toc
xabsr = @(x) rectanglefun(x, 75, 80) + rectanglefun(x, -80, -75)
figure
plot(x, xabsr(x))
[pdist, pint] = pdistW(U(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdist)
[pdista, pint] = pdistW(Ua(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdista)
hold on
plot(pint, pdista, 'r')
plot(pint, pdista*1e3, 'r')
plot(pint, pdista*1e2, 'r')
pi/dx
xabsr1 = @(x) rectanglefun(x, 79, 80) + rectanglefun(x, -80, -79)
[pdist1, pint] = pdistW(U(:, 1:20:end), xdomain, xabsr1);
figure
plot(pint, pdist1)
xabsr1 = @(x) rectanglefun(x, 75, 80).*(x-75).^2 + rectanglefun(x, -80, -75).*(x+75).^2
figure
plot(x, xabsr1(x))
[pdist1, pint] = pdistW(U(:, 1:20:end), xdomain, xabsr1);
figure
plot(pint, pdist1)
[U1 mniter matvecs] = TDHxpKr(K, Vf(x), @(u,x,t) -0.1*x*sin(0.01*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U1 mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(0.01*t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5);
viewWigner(U1(:,1:20:end), xdomain, 0.01, Wrangea
Wrange1 = getWrange(U1(:,1:20:end))
viewWigner(U1(:,1:20:end), xdomain, 0.01, Wrange1)
[pdist11, pint] = pdistW(U1(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdist11)
[U2 mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.1*x*sin(t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5);
mniter
viewVP(Ua(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
viewVP(U2(:,1:20:end), Vf, 0.1*sin((0:4:1e3)), x)
Wrange2 = getWrange(U2(:,1:20:end))
viewWigner(U2(:,1:20:end), xdomain, 0.01, Wrange2)
[pdist2, pint] = pdistW(U2(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdist2)
viewWigner(U(:,1:20:end), xdomain, 0.01, Wrange2)
[U3 mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.5*x*sin(0.01*t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5);
[U3 mniter matvecs] = TDHcheb_tsnf(Vf1, @(x,t) -0.2*x*sin(0.01*t), [-2, 15.5], fi0, xdomain, 1e3, 5e3, 9, 9, 1e-5);
mniter
viewVP(U3(:,1:20:end), Vf, 0.1*sin((0:4:1e3)), x)
viewVP(U3(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
viewVP(U3(:,1:20:end), Vf, 0.2*sin((0:4:1e3)*0.06), x)
viewVP(U3(:,1:20:end), Vf, 0.2*sin((0:4:1e3)*0.01), x)
Wrange3 = getWrange(U3(:,1:20:end))
viewWigner(U3(:,1:20:end), xdomain, 0.01, Wrange3)
[pdist3, pint] = pdistW(U3(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdist3)
viewPp(Ua(:,1:20:end), dx, xdlength, 0.01)
pabsr = @(x) rectanglefun(x, 75, 80).*(x-75).^2 + rectanglefun(x, -80, -75).*(x+75).^2
75/80*pi/dx
pabsr = @(x) rectanglefun(x, 4.71, pi/dx) + rectanglefun(x, -pi/dx, -4.71)
figure
plot(p, pabsr)
plot(p, pabsr(p))
[xdist, xint] = xdistW(U, xdomain, pabsr);
[xdist, xint] = xdistW(U(:, 1:20:end), xdomain, pabsr);
figure
plot(xint, xdist)
Vf
Kf = @(p) p.^2/2 + 1i*log(0.5*(tanh((p+4.71)*80*dx/pi)-tanh((p-75)*80*dx/pi)))
figure
plot(p, imag(Kf(p)))
p
Kf = @(p) p.^2/2 + 1i*log(0.5*(tanh((p+4.71)*80*dx/pi)-tanh((p-4.71)*80*dx/pi)))
plot(p, imag(Kf(p)))
[Uxp mniter matvecs] = TDHxpKr(Kf(p), Vf(x), @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
viewPp(Uxp(:,1:20:end), dx, xdlength, 0.01)
viewWigner(Uxp(:,1:20:end), xdomain, 0.01, Wrange
Wrangexp = getWrange(Uxp(:,1:20:end))
viewWigner(Uxp(:,1:20:end), xdomain, 0.01, Wrangexp)
figure
plot(pint, pdist1)
plot(pint, pdist)
plot(pint, pdist11)
plot(pint, pdista)
[pdistxp, pint] = pdistW(Uxp(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdistxp)
[xdista, xint] = pdistW(Ua(:, 1:20:end), xdomain, pabsr);
figure
plot(xint, xdista)
[xdistxp, xint] = pdistW(Uxp(:, 1:20:end), xdomain, pabsr);
figure
plot(xint, xdistxp)
max(abs(Ua(:,end)-Uxp(:,end)))
[Up mniter matvecs] = TDHxpKr(Kf(p), Vf1(x), @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
viewPp(Up(:,1:20:end), dx, xdlength, 0.01)
viewVP(Up(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
Kp
Kf
figure
plot(p, imag(Kf(p)))
figure
plot(x, imag(Vf(x)))
whos
%-- 29/12/2013 16:15 --%
load
whos
figure
Egrid = 0:188;
plot(Egrid, log(0.5*(1-tanh(x-50))))
plot(Egrid, log(0.5*(1-tanh(Egrid-50))))
plot(Egrid, log(0.5*(1-tanh(Egrid-50)+ 1-tanh(5.1))))
hold on
plot(Egrid, log(0.5*(1-tanh(Egrid-50))), 'r')
plot(Egrid, log(0.5*(1-tanh(Egrid-50)+ (1-tanh(5.1))*(tanh(51)+1))), 'g')
plot(Egrid, log(0.5*(1-tanh(Egrid-50)+ (1-tanh(5.1))*0.5*(tanh(Egrid-51)+1))), 'g')
plot(Egrid, log(0.5*(1-tanh(Egrid-50)+ (1-tanh(5.1))*0.5*(tanh(Egrid-55)+1))), 'g')
fabsE = @(x)  log(0.5*(1-tanh(x-50)+ (1-tanh(5.1))*0.5*(tanh(x-55)+1)))
whos
Ha = H + 1i*fabsE(H);
[Pa, Da] = eig(Ha);
Ea=diag(Da);
[Eareal, orderEa] = sort(real(Ea));
Ea = Ea(orderEa);
Pa = Pa(:, orderEa);
Ea0=Ea(1);
Ea(1:10)
Ea(1:50)
E(1:50)
[UE mniter matvecs] = TDHmat(Ha, @(u,t) XM*cos(t), [], [-1 188], fi0, [0 1e3], 5e3, 9, 9, 1e-5);
save
%-- 30/12/2013 12:53 --%
testTLSmat
profile
profile viewer
testTLSmat
load
clear all
testTLSmat
[U1 mniter matvecs] = TDHmatKr(H0, Vtfun, [], ui, [0 T], Nt, Nt_ts, Ncheb, tol);
mniter
matvecs
error1 = conj(U1(2, :)).*U1(2, :) - sin(t).^2;
maxer1 = max(abs(error1))
clear all
load
[UE mniter matvecs] = TDHmatKr(Ha, @(u,t) XM*cos(t), [], fi0, [0 1e3], 5e3, 9, 9, 1e-5);
[UE mniter matvecs] = TDHmatKr(Ha, @(u,t) XM*cos(t), [], fi0, [0 1e3], 1e4, 9, 9, 1e-5);
[UE mniter matvecs] = TDHmatKr(Ha, @(u,t) XM*cos(t), [], fi0, [0 1e3], 2e4, 9, 9, 1e-5);
[UE mniter matvecs] = TDHmatKr(Ha, @(u,t) XM*cos(t), [], fi0, [0 30], 600, 9, 9, 1e-5);
mniter
matvecs
viewVP(UE, Vf, -cos(0:0.05:30), x)
L/2
ans^2/2
E
fabsE = @(x)  log(0.5*(1-tanh(x-90)+ (1-tanh(5.1))*0.5*(tanh(x-95)+1)))
Ha = H + 1i*fabsE(H);
[UE mniter matvecs] = TDHmatKr(Ha, @(u,t) XM*cos(t), [], fi0, [0 30], 600, 9, 9, 1e-5);
mniter
viewVP(UE, Vf, -cos(0:0.05:30), x, 0.05)
viewWigner(UE, xdomain, 0.01, WrangeE)
WrangeE = getWrange(UE)
viewWigner(UE, xdomain, 0.01, WrangeE)
sqrt(180)
viewWigner(UE, xdomain, 0.01, WrangeE)
save
whos
viewWigner(UE, xdomain, 0.01, WrangeE)
viewWigner(Uxp, xdomain, 0.01, Wrangexp)
clear all
whos
Vf
Vf1 = @(x)1-1./sqrt(x.^2+1)
Kf = @(p) p.^2/2 + 1i*log(0.5*(tanh((p+4.71)*80*dx/pi)-tanh((p-4.71)*80*dx/pi)))
[Ux mniter matvecs] = TDHxpKr(K, Vf(x), @(u,x,t) -0.1*x*sin(0.01*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
min_x = xdomain(1);
max_x = xdomain(2);
xdlength = max_x - min_x;
dx = xdlength/Nx;
x = (min_x:dx:(max_x - dx)).';
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
[Ux mniter matvecs] = TDHxpKr(K, Vf(x), @(u,x,t) -0.1*x*sin(0.01*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
pabsr = @(x) rectanglefun(x, 4.71, pi/dx) + rectanglefun(x, -pi/dx, -4.71)
xabsr = @(x) rectanglefun(x, 75, 80) + rectanglefun(x, -80, -75)
[pdistx, pint] = pdistW(Ux(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdistx)
[Uxp mniter matvecs] = TDHxpKr(Kf(p), Vf(x), @(u,x,t) -0.1*x*sin(0.01*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
dx = xdlength/Nx;
Nx
[Uxp mniter matvecs] = TDHxpKr(Kf(p), Vf(x), @(u,x,t) -0.1*x*sin(0.01*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
dx
Kf = @(p) p.^2/2 + 1i*log(0.5*(tanh((p+4.71)*80*dx/pi)-tanh((p-4.71)*80*dx/pi)))
[Uxp mniter matvecs] = TDHxpKr(Kf(p), Vf(x), @(u,x,t) -0.1*x*sin(0.01*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
Wrangexp = getWrange(Uxp(:,1:20:end))
[pdistxp, pint] = pdistW(Uxp(:, 1:20:end), xdomain, xabsr);
figure
plot(pint, pdistxp)
hold on
viewVP(Ux(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
viewVP(Ux(:,1:20:end), Vf, 0.1*sin((0:4:1e3)), x)
viewVP(Ux(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.01), x)
hold on
plot(pint, pdistxp, 'r')
%-- 02/01/2014 14:54 --%
load
whos
clear all
load shlishi
whos
160*dx
dx
dx*256
clear all
load coulomb_abs80
whos
dx=160/256;
256-80
176/2
pi/dx
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', tol, 'TolX', tol);
tol=1e-10;
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', tol, 'TolX', tol);
whos
clear all
dx=160/256;
tol=1e-10;
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', tol, 'TolX', tol);
[
]
Vgx = [(rand(25,1) - 0.5); -rand(25,1)]
Vgx = [2*(rand(25,1) - 0.5); -rand(25,1)]
Vgx = [2*(rand(40,1) - 0.5); -rand(40,1)]
Vgx = [2*(rand(10,1) - 0.5); -rand(10,1)]
Vg = Vx2Vk0lb(Vgx)
[Vopt, optval, flag, data] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 25, 1e-10, 1, 1e3), Vg, options)
[Vopt1, optval1, flag1, data1] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 49, 1e-10, 1, 1e3), Vg, options)
[Vopt1, optval1, flag1, data1] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 49, 1e-10, 1, 1e3), Vopt, options)
figure
Vopt1x = dctI(Vopt1(1:11)) + 1i*dctI(Vopt1(12:22));
size(Vopt1)
Vopt1x = dctI0lb(Vopt1(1:10)) + 1i*dctI0lb(Vopt1(11:20));
plot(0.2:0.2:5, real(Vopt1x))
plot(0:dx:40*dx, real(Vopt1x))
plot(0:4*dx:40*dx, real(Vopt1x))
plot(0:4*dx:40*dx, imag(Vopt1x))
[Vopt2, optval2, flag2, data2] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 49, 1e-10, 1, 1e3), Vopt, options
[Vg2, Vg2x] = newguess0lb(Vopt1, 20);
figure
plot(0:2*dx:40*dx, imag(Vg2x))
plot(0:2*dx:40*dx, real(Vg2x))
[Vopt2, optval2, flag2, data2] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 49, 1e-10, 1, 1e3), Vg2, options)
[Vg3, Vg3x] = newguess0lb(Vopt2, 40);
figure
plot(0:dx:40*dx, real(Vg3x))
plot(0:dx:40*dx, imag(Vg3x))
plot(0:dx:40*dx, real(Vg3x))
[Vopt2, optval2, flag2, data2] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 49, 0.5e-10, 1, 1e3), Vopt2, options)
Vopt2x = dctI0lb(Vopt2(1:20)) + 1i*dctI0lb(Vopt2(21:40));
figure
plot(0:2dx:40*dx, real(Vg3x))
plot(0:2*dx:40*dx, real(Vopt2x))
hold on
plot(0:2*dx:40*dx, imag(Vopt2x), 'r')
40*dx-10
ans/dx
[Vopt2int, Vopt2xint] = newguess0lb(Vopt2, 40);
figure
plot(0:dx:40*dx, real(Vopt2xint))
plot(0:dx:40*dx, imag(Vopt2xint))
[Vopt2int, Vopt2xint] = newguess0lb_sp(Vopt2, 40);
plot(0:dx:40*dx, real(Vopt2xint))
figure
[allT, allR, allnorm] = allTRcoef1(Vopt2x,  [0 40*dx], [0.2 5], 97);
plot(0.2:0.05:5, allnorm)
[allT1, allR1, allnorm1] = allTRcoef1(Vopt1x,  [0 40*dx], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm1)
figure
plot(0:4*dx:40*dx, imag(Vopt1x))
plot(0:4*dx:40*dx, imag(Vg2x))
plot(0:2*dx:40*dx, imag(Vg2x))
plot(0:2*dx:40*dx, real(Vg2x))
[Vopt1int, Vopt1xint] = newguess0lb(Vopt1, 40);
figure
plot(0:dx:40*dx, real(Vopt1xint))
plot(0:dx:40*dx, imag(Vopt1xint))
optval
optval1
whos
[Voptp, optvalp, flagp, datap] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 25, 1e-12, 1, 1e3), Vg, options)
figure
max(abs(Vopt-Voptp))
plot(Vopt)
hold on
plot(Voptp, 'r')
[Voptp, optvalp, flagp, datap] = fminunc(@(V) percosV0lb_grad(V, [0 40*dx], [0.2 5], 25, 0, 1, 1e3), Vg, options)
plot(Voptp, 'g')
max(abs(Vopt-Voptp))
clear Voptp optvalp flagp datap
[Vopt24, optval24, flag24, data24] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg, options)
24*dx
Vgx24 = [2*(rand(12,1) - 0.5); -rand(12,1)]
Vg24 = Vx2Vk0lb(Vgx24)
[Vopt24, optval24, flag24, data24] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg24, options)
figure
Vopt24x = dctI0lb(Vopt1(1:12)) + 1i*dctI0lb(Vopt1(13:24));
Vopt24x = dctI0lb(Vopt24(1:12)) + 1i*dctI0lb(Vopt24(13:24));
plot(0:2*dx:24*dx, real(Vopt24x))
plot(0:2*dx:24*dx, imag(Vopt24x))
[Vopt24, optval24, flag24, data24] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vopt24, options)
Vopt24x = dctI0lb(Vopt24(1:12)) + 1i*dctI0lb(Vopt24(13:24));
figure
plot(0:2*dx:24*dx, imag(Vopt24x))
plot(0:2*dx:24*dx, real(Vopt24x))
plot(0:2*dx:24*dx, imag(Vopt24x))
Vgx24 = [2*(rand(12,1) - 0.5); -rand(12,1)]
[Vopt24, optval24, flag24, data24] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg24, options)
Vg24 = Vx2Vk0lb(Vgx24)
[Vopt24, optval24, flag24, data24] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg24, options)
Vgx24 = [2*(rand(12,1) - 0.5); -rand(12,1)]
Vg24 = Vx2Vk0lb(Vgx24);
[Vopt24, optval24, flag24, data24] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg24, options)
Vopt24x = dctI0lb(Vopt24(1:12)) + 1i*dctI0lb(Vopt24(13:24));
figure
plot(0:2*dx:24*dx, imag(Vopt24x))
plot(0:2*dx:24*dx, real(Vopt24x))
Vgx242 = [2*(rand(6,1) - 0.5); -rand(6,1)]
[Vopt242, optval242, flag242, data242] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg242, options)
Vg242 = Vx2Vk0lb(Vgx242);
[Vopt242, optval242, flag242, data242] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg242, options)
figure
Vopt242x = dctI0lb(Vopt242(1:6)) + 1i*dctI0lb(Vopt242(7:12));
plot(0:4*dx:24*dx, real(Vopt24x))
plot(0:4*dx:24*dx, real(Vopt242x))
plot(0:4*dx:24*dx, imag(Vopt242x))
Vgx242 = [2*(rand(6,1) - 0.5); -rand(6,1)];
Vg242 = Vx2Vk0lb(Vgx242);
[Vopt242, optval242, flag242, data242] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vg242, options)
[Vopt242, optval242, flag242, data242] = fminunc(@(V) percosV0lb_grad(V, [0 24*dx], [0.2 5], 25, 0, 1, 1e3), Vopt242, options)
figure
Vopt242x = dctI0lb(Vopt242(1:6)) + 1i*dctI0lb(Vopt242(7:12));
plot(0:4*dx:24*dx, imag(Vopt242x))
plot(0:4*dx:24*dx, real(Vopt242x))
32*dx
Vgx30 = [2*(rand(8,1) - 0.5); -rand(8,1)];
Vg30 = Vx2Vk0lb(Vgx242);
[Vopt30, optval30, flag30, data30] = fminunc(@(V) percosV0lb_grad(V, [0 30*dx], [0.2 5], 25, 0, 1, 1e3), Vopt30, options)
[Vopt30, optval30, flag30, data30] = fminunc(@(V) percosV0lb_grad(V, [0 30*dx], [0.2 5], 25, 0, 1, 1e3), Vg30, options)
Vgx30 = [2*(rand(8,1) - 0.5); -rand(8,1)];
Vg30 = Vx2Vk0lb(Vgx242);
Vg30 = Vx2Vk0lb(Vgx30);
[Vopt30, optval30, flag30, data30] = fminunc(@(V) percosV0lb_grad(V, [0 32*dx], [0.2 5], 25, 0, 1, 1e3), Vg32, options)
Vgx32 = [2*(rand(8,1) - 0.5); -rand(8,1)];
Vg32 = Vx2Vk0lb(Vgx32);
clear Vgx30 Vg30
[Vopt32, optval32, flag32, data32] = fminunc(@(V) percosV0lb_grad(V, [0 32*dx], [0.2 5], 25, 0, 1, 1e3), Vg32, options)
[Vopt32, optval32, flag32, data32] = fminunc(@(V) percosV0lb_grad(V, [0 32*dx], [0.2 5], 25, 0, 1, 1e3), Vopt32, options)
Vopt32x = dctI0lb(Vopt32(1:8)) + 1i*dctI0lb(Vopt32(9:16));
figure
plot(0:4*dx:30*dx, real(Vopt32x))
plot(0:4*dx:32*dx, real(Vopt32x))
plot(0:4*dx:32*dx, imag(Vopt32x))
4*dx
Vgx64 = [2*(rand(16,1) - 0.5); -rand(16,1)];
Vg64 = Vx2Vk0lb(Vgx32);
[Vopt64, optval64, flag64, data64] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 25, 0, 1, 1e3), Vopt64, options)
[Vopt64, optval64, flag64, data64] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 25, 0, 1, 1e3), Vg64, options)
[Vopt64, optval64, flag64, data64] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 25, 0, 1, 1e3), Vopt64, options)
figure
plot(0:4*dx:64*dx, imag(Vopt64x))
Vopt64x = dctI0lb(Vopt64(1:16)) + 1i*dctI0lb(Vopt64(17:32));
size(Vopt64)
Vg64 = Vx2Vk0lb(Vgx64);
[Vopt64, optval64, flag64, data64] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 25, 0, 1, 1e3), Vg64, options)
figure
plot(0:8*dx:64*dx, imag(Vopt64x))
Vopt64x = dctI0lb(Vopt64(1:8)) + 1i*dctI0lb(Vopt64(9:16));
plot(0:8*dx:64*dx, imag(Vopt64x))
plot(0:8*dx:64*dx, real(Vopt64x))
[Vg641, Vg641] = newguess0lb(Vopt64, 16);
[Vg641, Vg641x] = newguess0lb(Vopt64, 16);
figur
figure
plot(0:4*dx:64*dx, imag(Vg641x))
plot(0:4*dx:64*dx, real(Vg641x))
[Vopt641, optval641, flag641, data641] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 25, 0, 1, 1e3), Vg641, options)
[allT1, allR1, allnorm1] = allTRcoef1(Vopt64,  [0 64*dx], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm1)
[allT1, allR1, allnorm1] = allTRcoef1(Vopt64x,  [0 64*dx], [0.2 5], 97);
plot(0.2:0.05:5, allnorm1)
[allT2, allR2, allnorm2] = allTRcoef1(Vg641x,  [0 64*dx], [0.2 5], 97);
plot(0.2:0.05:5, allnorm2)
Vg641
percosV0lb_grad(Vg641, [0 64*dx], [0.2 5], 25, 0, 1, 1e3)
exp(Vg641x(19:end))
exp(imag(Vg641x(2:end)))
imag(Vg641x(2:end))
Vg641(2)
Vg641x(2)
Vg641x(2) = conj(Vg641x(2))
Vg641 = Vx2Vk0lb([real(Vg641x(2:end)); imag(Vg641(2:end))])
[real(Vg641x(2:end)); imag(Vg641(2:end))]
[real(Vg641x(2:end)); imag(Vg641x(2:end))]
Vg641 = Vx2Vk0lb([real(Vg641x(2:end)); imag(Vg641x(2:end))])
[Vopt641, optval641, flag641, data641] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 25, 0, 1, 1e3), Vg641, options)
Vopt641x = dctI0lb(Vopt641(1:16)) + 1i*dctI0lb(Vopt641(17:32));
[allT, allR, allnorm] = allTRcoef1(Vopt641x,  [0 64*dx], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm)
[Vopt642, optval642, flag642, data642] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 49, 0, 1, 1e3), Vg642, options)
[Vopt642, optval642, flag642, data642] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 49, 0, 1, 1e3), Vopt641, options)
save thursday
%-- 06/01/2014 12:15 --%
load thursday
whos
optval242
optval642
load coulomb_abs80
figure
plot(x, Vf(x))
size(Vopt642)
[Vint642, Vint642x] = newguess0lb(Vopt642, 64);
size
size(Vint642)
size(Vint642x)
figure
plot(0:dx:64*dx, imag(Vint642x))
plot(0:dx:64*dx, real(Vint642x))
[allT, allR, allnorm] = allTRcoef1(Vint642x,  [0 64*dx], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm)
plot(allR, allnorm)
plot(0.2:0.05:5, allR)
plot(0.2:0.05:5, allT)
[Vopt643, optval643, flag643, data643] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 97, 0, 1, 1e3), Vopt642, options)
[Vint643, Vint643x] = newguess0lb(Vopt642, 64);
figure
plot(0:dx:64*dx, real(Vint643x))
plot(0:dx:64*dx, imag(Vint643x))
[allT1, allR1, allnorm1] = allTRcoef1(Vint643x,  [0 64*dx], [0.2 5], 193);
figure
plot(0.2:0.025:5, allnorm1)
size(Vg64)
[Vopt644, optval644, flag644, data644] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 97, 0, 1, 1e3), Vg64, options)
[Vopt644, optval644, flag644, data644] = fminunc(@(V) percosV0lb_grad(V, [0 64*dx], [0.2 5], 97, 0, 1, 1e3), Vopt644, options)
optval643
optval644
save thursday
%-- 08/01/2014 14:59 --%
load thursday
mu = bound0der(@(x) x, x, -70, 70, 1);
figure
plot(x, mu)
%-- 09/01/2014 18:10 --%
integral(@(x) x^2, 0, 1)
integral(@(x) x.^2, 0, 1)
sech(5)
integral(@(x) x.^2, 1, 0)
mu = bound0der1(@(x) x, x, -70, 70, 1);
load thursday
mu = bound0der1(@(x) x, x, -70, 70, 1);
figure
plot(x, mu)
%-- 12/01/2014 18:22 --%
load thursday
whos
[Vint644, Vint644x] = newguess0lb(Vopt644, 64);
size(Vint644x)
size(Vint644)
figure
plot(0:dx:64*dx, real(Vint643x))
hold on
plot(0:dx:64*dx, real(Vint644x), 'r')
[allT2, allR2, allnorm2] = allTRcoef1(Vint644x,  [0 64*dx], [0.2 5], 193);
figure
plot(0.2:0.025:5, allnorm1)
plot(0.2:0.025:5, allnorm2)
hold on
plot(0.2:0.025:5, allnorm1, 'r')
optval644
optval643
Vabs_efficiency1(Vint643x, [0 64*dx], [0.2 5], 193)
Vabs_efficiency1(Vint644x, [0 64*dx], [0.2 5], 193)
Vabs_efficiency1(Vint644x, [0 64*dx], [0.2 5], 97)
Vabs_efficiency1(Vint643x, [0 64*dx], [0.2 5], 97)
[Vint644, Vint644x] = newguess0lb_sp(Vopt644, 64);
Vabs_efficiency1(Vint644x, [0 64*dx], [0.2 5], 97)
[Vint644, Vint644x] = newguess0lb(Vopt644, 64);
Vabs_efficiency1(Vint644x, [0 64*dx], [0.2 5], 97)
[Vintsp644, Vintsp644x] = newguess0lb_sp(Vopt644, 64);
[Vintsp643, Vintsp643x] = newguess0lb_sp(Vopt643, 64);
Vabs_efficiency1(Vint643x, [0 64*dx], [0.2 5], 97)
Vabs_efficiency1(Vintsp643x, [0 64*dx], [0.2 5], 97)
[allT2, allR2, allnorm2] = allTRcoef1(Vintsp644x,  [0 64*dx], [0.2 5], 193);
plot(0:dx:64*dx, real(Vintsp644x), 'g')
plot(0:dx:64*dx, real(Vintsp643x), 'm')
plot(0:dx:64*dx, real(Vint643x), 'm')
[Vint643, Vint643x] = newguess0lb(Vopt643, 64);
plot(0:dx:64*dx, real(Vint643x))
plot(0:dx:64*dx, real(Vintsp643x), 'm')
[allT1, allR1, allnorm1] = allTRcoef1(Vint643x,  [0 64*dx], [0.2 5], 193);
figure
Vabs_efficiency1(Vintsp643x, [0 64*dx], [0.2 5], 97)
Vabs_efficiency1(Vint643x, [0 64*dx], [0.2 5], 97)
plot(0.2:0.025:5, allnorm1, 'r')
[allT1, allR1, allnorm1] = allTRcoef1(Vint643x,  [0 64*dx], [0.2 5], 193);
plot(0.2:0.025:5, allnorm1, 'r')
[allT3, allR3, allnorm3] = allTRcoef1(Vintsp643x,  [0 64*dx], [0.2 5], 193);
plot(0.2:0.025:5, allnorm3, 'g')
Vabs_efficiency1(Vintsp643x, [0 64*dx], [0.2 5], 97)
Vabs_efficiency1(Vintsp644x, [0 64*dx], [0.2 5], 97)
Vabs_efficiency1(Vintsp644x, [0 64*dx], [0.2 5], 193)
Vabs_efficiency1(Vintsp643x, [0 64*dx], [0.2 5], 193)
figure
plot(0:dx:64*dx, imag(Vintsp643x))
hold on
plot(0:dx:64*dx, imag(Vintsp644x), 'r')
save thursday
%-- 13/01/2014 15:54 --%
%-- 15/01/2014 12:04 --%
load thursday
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3), Vg64, options
whos
Vg64
Vg64x
Vgx64
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', tol, 'TolX', tol, 'MaxIter', 1e4);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e10), Vg64, options)
Vgx645 = [2*(rand(16,1) - 0.5); -rand(16,1)-1];
Vg645 = Vx2Vk0lb(Vgx645)
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e10), Vg645, options);
Vgx645 = [2*(rand(16,1) - 0.5); -rand(16,1)];
Vg645 = Vx2Vk0lb(Vgx645);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e10), Vg645, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3), Vg645, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 49, 0, 1, 1e10), Vg645, options
pi/dx
0.2/(2*pi)*dx
2*pi/xdlength
2*pi/160
min_x = xdomain(1);
max_x = xdomain(2);
xdlength = max_x - min_x;
dx = xdlength/Nx;
x = (min_x:dx:(max_x - dx)).';
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
size(p)
p(1:128)
dp = 2*pi/xdlength;
p(129)
p(7)
p(6)
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [p(6) p(128)], 64, 49, 0, 1, 1e10), Vg645, options
size(p(6:2:128))
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [p(6) p(128)], 64, 62, 0, 1, 1e10), Vg645, options);
Vgx645 = [2*(rand(16,1) - 0.5); -rand(16,1)];
Vg645 = Vx2Vk0lb(Vgx645);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [p(6) p(128)], 64, 62, 0, 1, 1e10), Vg645, options);
Vgx645 = [2*(rand(16,1) - 0.5); -rand(16,1)];
Vg645 = Vx2Vk0lb(Vgx645);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [p(6) p(128)], 64, 62, 0, 1, 1e10), Vg645, options);
Vgx645 = [2*(rand(16,1) - 0.5); -rand(16,1)*10];
Vg645 = Vx2Vk0lb(Vgx645);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [p(6) p(128)], 64, 62, 0, 1, 1e10), Vg645, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e10), Vg645, options);
[performance, gradient] = percosV0lb_grad1(Vopt644, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3);
performance
[performance, gradient] = percosV0lb_grad1(Vopt644, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3);
Vximag(0:4:end)
Vximag(1:4:end)
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e10), Vg645, options);
figure
plot(0:4*dx:64*dx, real(Vg644x))
plot(0:4*dx:64*dx, real(Vopt644x))
plot(0:4*dx:64*dx, real(Voptx644))
whos
plot(0:dx:64*dx, real(Vint644x))
[performance, gradient] = percosV0lb_grad1(Vopt644, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3);
hold on
plot(0:dx:64*dx, Vxreal, 'r')
dx=160/256
plot(0:dx:64*dx, Vxreal, 'r')
plot(0:dx:64*dx, Vxreal*2, 'r')
plot(0:dx:64*dx, Vximag*2, 'g')
gradnorms
grad_norms
grad_perf_Vx
figure
plot(0:dx:64*dx, Vximag)
hold on
plot(0:dx:64*dx, log(gradient(66:130)), 'r')
plot(0:dx:64*dx, log(grad_perf_Vx(66:130)), 'r')
plot(0:dx:64*dx, log(grad_perf_Vx(66:130))*1e-2, 'r')
figure
plot(0:dx:64*dx, Vxreal)
hold on
plot(0:dx:64*dx, grad_perf_Vx(1:65), 'r')
plot(0:dx:64*dx, grad_perf_Vx(1:65)*1e-1, 'r')
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e10), Vg645*2, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), Vg645*2, options);
flag645
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 49, 0, 1, 1e3), Vgopt645, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 49, 0, 1, 1e3), Vopt645, options);
%-- 16/01/2014 04:19 --%
%-- 16/01/2014 19:18 --%
load thursday
[performance, gradient] = percosV0lb_grad2(Vopt644, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3);
performance
[performance1, gradient1] = percosV0lb_grad1(Vopt644, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3);
performance1
sqrt(performance)
[performance1, gradient1] = percosV0lb_grad1(2*Vopt644, [0 64*dx], [0.2 5], 64, 97, 0, 1, 1e3);
performance1
[performance, gradient] = percosV0lb_grad2(Vopt644, [0 64*dx], [0.2 5], 16, 97, 0, 1, 1e3);
performance1
performance
Vopt644val
optval644
Vgx645 = [2*(rand(16,1) - 0.5); -rand(16,1)];
Vg645 = Vx2Vk0lb(Vgx645);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad2(V, [0 64*dx], [0.2 5], 64, 49, 0, 1, 1e3), Vopt645, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), Vg645, options);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad2(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), Vg645, options);
figure
plot(0:4*dx:64*dx, Vximag)
dx=160/256
plot(0:4*dx:64*dx, Vximag)
hold on
plot(0:dx:64*dx, Vximag_penal, 'r')
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', tol, 'TolX', tol, 'MaxIter', 1e4);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad2(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), Vopt645, options);
options = optimoptions(@fminunc,'GradObj','on','PlotFcns',{@optimplotfval,@optimplotfirstorderopt}, 'TolFun', tol, 'TolX', tol, 'MaxIter', 1e4, 'MaxFunEvals', 1e5);
[Vopt645, optval645, flag645, data645] = fminunc(@(V) percosV0lb_grad2(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), Vopt645, options);
Vopt645x = dctI0lb(Vopt645(1:16)) + 1i*dctI0lb(Vopt645(17:32));
Vabs_efficiency1(Vopt645x, [0 64*dx], [0.2 5], 25)
Vabs_efficiency1(Vopt645x, [0 64*dx], [0.2 5], 49)
[Vopt646, optval646, flag646, data646] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), 2*Vopt645, options);
[Vopt646, optval646, flag646, data646] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), Vopt645, options);
[Vopt646, optval646, flag646, data646] = fminunc(@(V) percosV0lb_grad1(V, [0 64*dx], [0.2 5], 64, 25, 0, 1, 1e3), 2*Vopt645, options);
optval645
optval646
save thursday
%-- 21/01/2014 12:37 --%
%-- 21/01/2014 12:39 --%
load thursday
whos
clear Vgx30 Vg30
clear all
load thursday
whos
Vf
Vf1
Vf1 = @(x)1-1./sqrt(x.^2+1)
figure
plot(x, Vf1(x))
256+128
x0d = bound0der1(@(x) x, x, -37.5, 37.5, 1);
figure
plot(x, x0d)
V0d = bound0der1(Vf1, x, -37.5, 37.5, 1);
figure
plot(x, V0d)
figure
plot(0:dx:64*dx, Vopt649)
optval649
[Vopt649, optval649, flag649, data649] = fminunc(@(V) percosV0lb_grad3(V, [0 64*dx], [0.2 5], 64, 49, 0), Vg649, options
size(Vopt648)
optval648
[Vopt649, optval649, flag649, data649] = fminunc(@(V) percosV0lb_grad3(V, [0 64*dx], [0.2 5], 64, 49, 0), Vopt648*2, options);
[~, Vopt649x] = newguess0lb3(Vopt649, 64);
figure
Vopt649x = Vopt649x/2;
plot(0:dx:64*dx, Vopt649x)
Vabs_efficiency1(Vopt649x, [0 64*dx], [0.2 5], 49)
Vopt649x = Vopt649x*2;
Vopt649x = Vopt649x/2;
Vabs_efficiency1(Vopt649x, [0 64*dx], [0.2 5], 49)
size(Vopt649
size(Vopt649)
size(Vopt649x)
run('C:\Users\Ido\Dropbox\MATLAB\interpct.m')
Vabs_efficiency1(Vopt649x*2, [0 64*dx], [0.2 5], 49)
Vabs_efficiency1(Vopt649x*4, [0 64*dx], [0.2 5], 49)
Vabs_efficiency1(Vopt649x/4, [0 64*dx], [0.2 5], 49)
Vabs_efficiency1(Vopt649x/2, [0 64*dx], [0.2 5], 49)
Vabs_efficiency1(Vopt649x, [0 64*dx], [0.2 5], 49)
[Vopt649, optval649, flag649, data649] = fminunc(@(V) percosV0lb_grad3(V, [0 64*dx], [0.2 5], 64, 49, 0), Vopt649, options);
hold on
plot(0:dx:64*dx, V)
plot(0:160/256:64*160/256, V)
figure
plot(0:160/256:64*160/256, imag(V))
hold on
plot(0:dx:64*dx, imag(Vopt649x), 'r')
plot(0:dx:64*dx, imag(Vopt649x)/4, 'g')
clear Vopt649x
save thursday
%-- 22/01/2014 12:18 --%
load thursday
Vopt649x = Vopt2Vx0lb3(Vopt649, 64);
figure
plot(0:dx:64*dx, real(Vopt649x))
figure
plot(0:dx:64*dx, imag(Vopt649x))
Vabs_efficiency1(Vopt649x, [0 64*dx], [0.2 5], 49)
[allT1, allR1, allnorm1] = allTRcoef1(Vopt649x,  [0 64*dx], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm)
dctimVopt649x = dctI(imag(Vopt649x));
figure
plot(0:pi/40:pi/dx, dctimVopt649x)
pi/(2*dx)
save thursday
%-- 26/01/2014 14:50 --%
load thursday
whos
optval6412
Vopt6412x = Vopt2Vx0lb3(Vopt6412, 64);
figure
plot(0:dx:64*dx, imag(Vopt6412x))
plot(0:dx:64*dx, real(Vopt6412x))
plot(0:dx:64*dx, imag(Vopt6412x))
dctimVopt6412x = dctI(imag(Vopt6412x));
figure
plot(0:pi/40:pi/dx, dctimVopt6412x)
pi/(2*dx)
Vabs_efficiency1(Vopt6412x, [0 64*dx], [0.2 5], 49)
Vabs_efficiency1(Vopt6412x, [0 64*dx], [0.2 5], 25)
Vabs_efficiency1(Vopt6412x, [0 64*dx], [0.2 5], 49)
Vopt6412
save thursday
%-- 26/01/2014 18:30 --%
load thursday
Voptval6412
optval6412
Vabs_efficiency1(Vopt6411x, [0 64*dx], [0.2 5], 49)
Vabs_efficiency1(Vopt6410x, [0 64*dx], [0.2 5], 49)
Vopt6412x = Vopt2Vx0lb4(Vopt6412, 64);
Vabs_efficiency1(Vopt6412x, [0 64*dx], [0.2 5], 49)
figure
plot(0:dx:64*dx, real(Vopt6412x))
plot(0:dx:64*dx, imag(Vopt6412x))
dctimVopt6412x = dctI(imag(Vopt6412x));
figure
plot(0:pi/40:pi/dx, dctimVopt6412x)
save thursday
[allT2, allR2, allnorm2] = allTRcoef1(Vopt6412x,  [0 64*dx], [0.2 5], 97);
figure
plot(0.2:0.05:5, allnorm)
plot(0.2:0.05:5, allnorm2)
[allT1, allR1, allnorm1] = allTRcoef1(Vopt649x,  [0 64*dx], [0.2 5], 97);
hold on
plot(0.2:0.05:5, allnorm1, 'r')
optval649
save thursday
%-- 27/01/2014 11:40 --%
load thursday
optval648
optval647
optval646
optval644
save optVabs Vopt649x Vopt6412x
optval649
optval6412
clear all
load optVabs
whos
Vf = @(x)1-1./sqrt(x.^2+1)
[fi0, E0, x, E, P, H] = gsV(Vf, xdomain, Nx);
xdomain = [-40 40];
Nx=128;
min_x = xdomain(1);
max_x = xdomain(2);
xdlength = max_x - min_x;
dx = xdlength/Nx;
x = (min_x:dx:(max_x - dx)).';
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
dx
K = p.^2/2;
xdomain = [-80 80];
Nx=256;
min_x = xdomain(1);
max_x = xdomain(2);
xdlength = max_x - min_x;
dx = xdlength/Nx;
x = (min_x:dx:(max_x - dx)).';
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
xabs = bound0der1(@(x) x, x, -37.5, 37.5, 1);
figure
plot(x, xabs)
V0 = bound0der1(Vf, x, -37.5, 37.5, 1);
figure
plot(x, V0)
256-65
Vabs = [V0(1:65)+Vopt649x(65:-1:1); V0(66:191); V0(192:256) + Vopt649x];
size(Vabs)
figure
plot(x, real(Vabs))
plot(x, imag(Vabs))
Vabs1 = [V0(1:65)+Vopt6412x(65:-1:1); V0(66:191); V0(192:256) + Vopt6412x];
plot(x, real(Vabs1))
plot(x, imag(Vabs1))
[fi0, E0, x, E, P, H] = gsV(Vabs, xdomain, Nx);
E(1:20)
[fi01, E01, x1, E1, P1, H1] = gsV(Vabs1, xdomain, Nx);
[E E1]
max(imag(Vabs1))
save Vabs_comparison
figure
plot(x, P(:,1))
max(abs(imag(P(:,1))))
max(abs(imag(P(:,2))))
max(abs(imag(P(:,3))))
max(abs(imag(P(:,30))))
max(abs(imag(P(:,200))))
max(abs(fi0-fi01))
max(abs(P(:,2)-P1(:,2)))
max(abs(real(P(:,2)-P1(:,2))))
max(abs(real(P(:,1)-P1(:,1))))
plot(x, P(:,2))
plot(x, P(:,3))
plot(x, P(:,4))
figure
plot(x, P1(:,3))
plot(x, P1(:,4))
plot(x, imag(P1(:,4)))
plot(x, P(:,5))
plot(x, P(:,6))
plot(x, P(:,7))
plot(x, P(:,8))
plot(x, P(:,9))
plot(x, P(:,10))
E(1:20)
[(1:20)', E(1:20)]
plot(x, P(:,12))
plot(x, P(:,13))
E1(1:20)
[(1:20)', E1(1:20)]
plot(x, P1(:,5))
plot(x, P1(:,6))
plot(x, P1(:,7))
plot(x, P1(:,8))
plot(x, P1(:,9))
size(Vopt6412x)
plot(x, P(:,100))
plot(x, P(:,130))
plot(x, P(:,200))
figure
plot(x, Vabs1)
size(V0)
size(V0abs)
size(Vabs)
Vabs = [V0(1:65)+Vopt649x(65:-1:1); V0(66:192); V0(193:256) + Vopt649x(1:64)];
figure
plot(x, Vabs1)
plot(x, Vabs)
Vabs1 = [V0(1:65)+Vopt6412x(65:-1:1); V0(66:192); V0(193:256) + Vopt6412x(1:64)];
[fi01, E01, x1, E1, P1, H1] = gsV(Vabs1, xdomain, Nx);
[fi0, E0, x, E, P, H] = gsV(Vabs, xdomain, Nx);
[E(1:20) E1(1:20)]
plot(x, P1(:,3))
plot(x, P1(:,4))
plot(x, P1(:,3))
plot(x, Vabs1)
plot(x, P1(:,4))
plot(x, conj(P1(:,4)).*P1(:,4))
plot(x, imag(P1(:,4)))
plot(x, imag(P1(:,5)))
plot(x, (P1(:,5)))
plot(x, (P1(:,6)))
plot(x, (P1(:,5)))
plot(x, (P1(:,6)))
plot(x, (P1(:,7)))
plot(x, (P1(:,6)))
plot(x, (P(:,13)))
[E, E1]
[(1:256)' E, E1]
figure
plot(x, (P(:,214)))
plot(x, (P(:,215)))
plot(x, (P1(:,145)))
plot(x, (P1(:,146)))
plot(x, (P1(:,145)))
save Vabs_comparison
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
viewVP(U(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
mniter
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs1, @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
vie
viewVP(U1(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
figure
plot(x, U(:, end))
hold on
plot(x, U1(:, end), 'r')
figure
plot(x, conj(U(:, end)).*U(:.end))
plot(x, conj(U(:, end)).*U(:,end))
hold on
plot(x, conj(U1(:, end)).*U1(:,end), 'r')
figure
plot(x, conj(U(:, end)).*U(:,end)-conj(U1(:, end)).*U1(:,end))
mniter1
mx=evx(U, x);
figure
plot(0:0.2:1e3, mx)
Kre = 250;
mx1=evx(U1, x);
c=mx(1:5:end);
size(c)
dt=1;
FBDMrinp
Kre = 125;
FBDMrinp
Kre = 250;
FBDMrinp1
c=mx(1:5:end).';
FBDMrinp
c=mx1(1:5:end).';
FBDMrinp
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
figure
plot(x, xabs)
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
U=U(:, 1:2:end);
viewVP(U(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), xabs)
viewVP(U(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
norm(U(:, end))
mniter
x1=x;
[U2 mniter2 matvecs2] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*x1*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
clear U U1
clear U U1 mx mx1 c
save Vabs_comparison
%-- 28/01/2014 12:26 --%
load Vabs_comparison
whos
figure
plot(x, x1)
plot(x, xabs)
Vtfun = @(u,x,t) -0.1*x1*sin(0.06*t);
figure
plot(x, Vtfun(1, 1, 0))
plot(x, Vtfun(1, 1, pi/0.12))
plot(x, Vtfun(1, 1, pi/0.24))
plot(x, Vtfun(1, 1, pi/0.03))
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
plot(x, Vtfun(1, 1, pi/0.12))
plot(x, Vtfun(1, x, pi/0.12))
plot(x, Vtfun(Ulast(:, tmidi), x, pi/0.12))
plot(x, Vtfun(Ulast(:, tmidi), x, t(tmidi)))
figure
plot(x, Vthalf)
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*x1*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
norm(U(:, end))
mniter
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
norm(U(:, end))
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*x*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
t
K
K = p.^2/2;
Kp = K;
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
viewVP(U(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs1, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
mniter1
viewVP(U1(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x)
figure
plot(x, U(:, end))
plot(x, conj(U(:, end)).*U(:,end))
hold on
plot(x, conj(U1(:, end)).*U1(:,end), 'r')
figure
plot(x, conj(U(:, end)).*U(:,end)-conj(U1(:, end)).*U1(:,end))
[wreal1, d_result1, w_result1, g_result1, meand1, w_maxd1, maxd1] = KBDMrinp(c, dt, Kre);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(c, dt, Kre);
mx=evx(U, x);
mx1=evx(U1, x);
c=mx(1:5:end);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(c, dt, Kre);
c=mx(1:5:end).';
w_maxd
FBDMrinp
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(mx, dt, Kre);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(c, dt, Kre);
[wreal1, d_result1, w_result1, g_result1, meand1, w_maxd1, maxd1] = KBDMrinp(c1, dt, Kre);
mx1=evx(U1, x);
size(U1)
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs1, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
K = p.^2/2;
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs1, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mx1=evx(U1, x);
c1=mx1(1:5:end).';
[wreal1, d_result1, w_result1, g_result1, meand1, w_maxd1, maxd1] = KBDMrinp(c1, dt, Kre);
figure
plot(wreal, meand)
plot(wreal(meand<1e3), meand(meand<1e3))
E(2)-E(1)
figure
plot(wreal1(meand1<1e3), meand1(meand1<1e3))
4.402/0.2781
4.021/0.2847
figure
plot(wreal, g_result)
plot(wreal, log(g_result))
xdomaine = [-160 160];
save Vabs_comparison
%-- 03/02/2014 22:28 --%
%-- 04/02/2014 14:12 --%
load Vabs_comparison
whos
clear B Du U0 U0cti
save Vabs_comparison
figure
plot(x, conj(Ue(:, end)).*Ue(:,end))
plot(xe, conj(Ue(:, end)).*Ue(:,end))
hold on
plot(xe2, conj(Ue2(:, end)).*Ue2(:,end), 'r')
xdomaine
xdomaine2
xdomaine3 = [-720 720];
0.06/(2*pi)
1/ans
5.33e-9*sqrt(2e-13)
5.33e-9*sqrt(2e13)
field13 = ans;
1/2.42e-17
ans/(2*pi)
3e8/ans
*1e9
ans*1e9
ans/(2*pi*137)
1/45.6
45.6/1064
1/(2*pi*137)
1/ans
1/5.29e-11
ans*86.1e-2
1.1617e-03*1e-2*1.8904e+10
1/ans
ans*2e3
0.06/4.55e-6
0.06/4.5536e-06
field13/0.429^2
field13
field13/0.0429^2
45.6/1064
%-- 05/02/2014 14:21 --%
whos
max(abs(fieldt))
(ans/5.33e-9)^2
clear all
max(abs(fieldt))
(ans/5.33e-9)^2
%-- 06/02/2014 15:13 --%
(0.1/5.33e-9)^2
load Vabs_comparison
whos
max(abs(xe21-xe2))
min(abs(xe21-xe2))
clear xe21
mxe2=evx(Ue2, xe2);
figure
plot(0:0.02:1e3, mxe2)
size(mxe2)
plot(0:0.2:1e3, mxe2)
0.1/0.06^2
ce2=mxe2(1:5:end).';
figure
plot(wreal, meande)
plot(wreale, meande)
xdomaine2
plot(wreale(wreale<100), meande(wreale<100))
plot(wreale(meande<100), meande(meande<100))
plot(wreale(meande<100), log(meande(meande<100)))
figure
plot(wreale, log(meande))
[wreale2, d_resulte2, w_resulte2, g_resulte2, meande2, w_maxde2, maxde2] = KBDMrinp(ce2, dt, Kre);
figure
plot(wreale2, meande2)
plot(wreale2(meande2<100), log(meande2(meande2<100)))
hold on
plot(wreale(meande<100), log(meande(meande<100)), 'r')
hold on
plot(wreale(meande<100), log(meande(meande<100)), 'r')
plot(wreal(meand<100), log(meand(meand<100)), 'g')
figure
plot(wreale(meande<100), (meande(meande<100)))
hold on
plot(wreale2(meande2<100), (meande2(meande2<100)), 'g')
xdomaine3 = [-640 640];
dx
xe3= (-640:dx:(640-dx)).';
size(xe3)
whos
Ve30 = bound0der1(Vf, x, -597.5, 597.5, 1);
xe3abs = bound0der1(@(x) x, x, -597.5, 597.5, 1);
figure
plot(xe3, xe3abs)
xe3abs = bound0der1(@(x) x, xe3, -597.5, 597.5, 1);
Ve30 = bound0der1(Vf, xe3, -597.5, 597.5, 1);
plot(xe3, xe3abs)
Ve3abs = [Ve30(1:65)+Vopt649x(65:-1:1); Ve30(66:1984); Ve30(1985:2048) + Vopt649x(1:64)];
figure
plot(x, Ve30)
plot(xe3, Ve30)
figure
plot(xe3, Ve3abs)
plot(xe3, imag(Ve3abs))
fi0e3 = [zeros(512,1); fi0e2; zeros(512, 1)];
figure
plot(xe3, fi0e3)
[Ue3 mnitere3 matvecse3] = TDHxpKr(Ke3, Ve3abs, @(u,x,t) -0.1*xe3abs*sin(0.06*t), [], fi0e3, xe3, [0 1e3], 2e4, 9, 9, 1e-5);
p = (0:(2*pi/1280):(2*pi*(1/dx - 1/1280))).';
xdlengthe3=2*640
Nxe
Nxe3=2048;
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
pe3 = (0:(2*pi/xdlengthe3):(2*pi*(1/dx - 1/xdlengthe3))).';
pe3((Nxe3/2 + 1):Nxe3) = pe3((Nxe3/2 + 1):Nxe3) - 2*pi/dx;
Ke3 = pe3.^2/2;
[Ue3 mnitere3 matvecse3] = TDHxpKr(Ke3, Ve3abs, @(u,x,t) -0.1*xe3abs*sin(0.06*t), [], fi0e3, xe3, [0 1e3], 2e4, 9, 9, 1e-5);
mnitere3
size(Ue3)
Ue3=Ue3(:, 1:4:end);
size(Ue3)
viewVP(Ue3(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), xe3)
viewPp(Ue3(:,1:20:end), dx, xdlengthe3, 0.01)
figure
plot(xe2, conj(Ue2(:, end)).*Ue2(:, end))
hold on
plot(xe3, conj(Ue3(:, end)).*Ue3(:, end))
plot(xe3, conj(Ue3(:, end)).*Ue3(:, end), 'r')
figure
norm(Ue2(:, end))
norm(Ue3(:, end))
plot(xe2, conj(Ue3(513:1554, end)).*Ue2(513:1554, end) - conj(Ue2(:, end)).*Ue2(:, end))
plot(xe2, conj(Ue3(513:1554, end)).*Ue3(513:1554, end) - conj(Ue2(:, end)).*Ue2(:, end))
2048-512
plot(xe2, conj(Ue3(513:1536, end)).*Ue3(513:1536, end) - conj(Ue2(:, end)).*Ue2(:, end))
mxe3=evx(Ue3, xe3);
ce3=mxe3(1:5:end).';
[wreale3, d_resulte3, w_resulte3, g_resulte3, meande3, w_maxde3, maxde3] = KBDMrinp(ce3, dt, Kre);
figure
plot(wreale3(meande3<100), (meande3(meande3<100)))
plot(wreale3(meande3<100), log(meande3(meande3<100)))
plot(wreale3(meande3<100), log(meande3(meande3<100)), 'm')
whos
xdomain240 = [-240 240];
xdlength240=2*240
x240= (-240:dx:(240-dx)).';
Nx240=xdlength240/dx
xdlengthe
Nxe
512+4*64
x240abs = bound0der1(@(x) x, x240, -197.5, 197.5, 1);
V0240 = bound0der1(Vf, x240, -197.5, 197.5, 1);
figure
plot(x240, x240abs)
plot(x240, V0240abs)
plot(x240, V0240)
V240abs = [Ve30(1:65)+Vopt649x(65:-1:1); V0240(66:704); V0240(705:768) + Vopt649x(1:64)];
figure
plot(x240, V240abs)
p240 = (0:(2*pi/xdlength240):(2*pi*(1/dx - 1/xdlength240))).';
p240((Nx240/2 + 1):Nx240) = p240((Nx240/2 + 1):Nx240) - 2*pi/dx;
[U240 mniter240 matvecs240] = TDHxpKr(K240, V240abs, @(u,x,t) -0.1*x240abs*sin(0.06*t), [], fi0240, x240, [0 1e3], 1e5, 9, 9, 1e-5);
K240 = p240.^2/2;
fi0240 = gsV(V240abs, xdomain240, Nx240);
figure
plot(x240, fi0240)
[U240 mniter240 matvecs240] = TDHxpKr(K240, V240abs, @(u,x,t) -0.1*x240abs*sin(0.06*t), [], fi0240, x240, [0 1e3], 1e5, 9, 9, 1e-5);
[U240 mniter240 matvecs240] = TDHxpKr(K240, V240abs, @(u,x,t) -0.1*x240abs*sin(0.06*t), [], fi0240, x240, [0 1e3], 5e3, 9, 9, 1e-5);
mniter240
norm(U240(:, end))
viewVP(U240(:,1:20:end), Vf, 0.1*sin((0:4:1e3)*0.06), x240)
figure
plot(xe3, fi0e3)
plot(xe3, conj(Ue3(:, end)).*Ue3(:, end))
hold on
plot(x240, conj(U240(:, end)).*U240(:, end), 'r')
figure
plot(x240, conj(Ue3(513:1536, end)).*Ue3(513:1536, end) - conj(Ue2(:, end)).*Ue2(:, end))
1024-768/2
1024+768/2
plot(x240, conj(Ue3(641:1408, end)).*Ue3(641:1408, end) - conj(U240(:, end)).*U240(:, end))
figure
plot(x240, (conj(U240(:, end)).*U240(:, end) - conj(Ue3(641:1408, end)).*Ue3(641:1408, end))./(conj(Ue3(641:1408, end)).*Ue3(641:1408, end)))
figure
plot(xe2, (conj(Ue2(:, end)).*Ue2(:, end) - conj(Ue3(513:1536, end)).*Ue3(513:1536, end))./(conj(Ue3(513:1536, end)).*Ue3(513:1536, end)))
whos
clear H H1 He He2 He21 Pe2 Pe21
save Vabs_comparison
%-- 12/02/2014 12:49 --%
whos
Vf = @(x)1-1./sqrt(x.^2+1)
xdomain = [-40 40];
Nx=128;
xdlength=80;
dx=xdlength/Nx
[fi0, E0, x, E, P, H] = gsV(Vf, xdomain, Nx);
E(1:10)
xdomain = [-40 40];
xdomain = [-80 80];
xdlenght=160;
load wed
whos
clear all
load thursday
whos
save optV40 Vopt649x
clear all
load optV40
whos
Vf = @(x)1-1./sqrt(x.^2+1)
xdomain = [-80 80];
xdlenght=160;
Nx=256;
dx=xdlength/Nx
xdlength=160;
clear xdlenght
dx=xdlength/Nx
whos
[Vabs, xabs, x, p, K, Nx] = get_prop_vars(Vf, xdomain, dx, Vopt649x);
[fi0, E0, x, E, P, H] = gsV(Vabs, xdomain, Nx);
E(1:10)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, miux, @(x) -10*(1 - x)^2, @(x) 20*(1 - x), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), orthpenal, 0, 1e-5, 2e3, 0.2, 7, 7, 1e-3);
figure
plot(x, xabs)
plot(x, Vabs)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) -10*(1 - x)^2, @(x) 20*(1 - x), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 2e3, 0.2, 7, 7, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) -10*(1 - x)^2, @(x) 20*(1 - x), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3);
mallniterc
norm(psi(:, end))
figure
plot(0:pi/1e3:pi/0.2, fieldw)
plot(0:pi/1e3:pi/1e1, fieldw(1:101))
viewVP(psi(:,1:20:end), Vf, fieldt, x)
viewVP(psi, Vf, fieldt, x)
viewVP(psi, Vf, fieldt, x, 0.01)
figure
plot(0:pi/1e3:pi/1e1, evmiuw(1:101))
plot(0:pi/1e3:pi/5, evmiuw(1:201))
plot(0:pi/1e3:pi, evmiuw(1:1001))
J1
conv(end)
figure
plot(0:0.2:1e3, fieldt)
2*pi/0.03
plot(0:0.2:1e3, evmiut)
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) -(1 - x)^2, @(x) 2*(1 - x), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 100);
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) -(1 - x)^2, @(x) 2*(1 - x), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(psi1(:, end))
norm(psi(:, end))
J1
figure
plot(0:0.2:1e3, evmiut)
plot(0:0.2:1e3, fieldt)
plot(0:0.2:1e3, (fieldt/5.33e-9).^2)
figure
plot(0:0.2:1e3, fieldt1)
J11
%-- 16/02/2014 14:17 --%
%-- 18/02/2014 10:17 --%
5.9
0.87
%-- 18/02/2014 11:36 --%
7.9
%-- 18/02/2014 11:41 --%
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
load('coulomb_optV40.mat')
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
67
9.8
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
norm(psi(:, end))
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-4, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 10);
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) -(1 - x)^2, @(x) 2*(1 - x), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 10);
norm(psi1(:, end))
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1);
viewVP(psi1, Vf, fieldt1, x, 0.01)
figure
plot(0:0.2:1e3, fieldt1)
plot(0:0.2:1e3, evmiut1)
plot(0:pi/1e3:pi, evmiuw(1:1001))
plot(0:pi/1e3:pi, evmiuw1(1:1001))
doc fminunc
doc optimplotfval
edit optimplotfval
edit fminunc
clear all
load('thursday.mat')
whos
[Vopt649, optval649, flag649, data649] = fminunc(@(V) percosV0lb_grad3(V, [0 64*dx], [0.2 5], 64, 49, 0), Vopt649, options);
doc setappdata
clear all
load('coulomb_optV40.mat')
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
norm(psi(:, end))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
figure
but = uicontrol('style', 'push', 'string', 'Stop');
get(but, 'position')
but1 = uicontrol('style', 'push', 'string', 'Stop1', 'position', [20, 0, 60, 20]);
set(but1, 'position', [20, -5, 60, 20])
set(but1, 'position', [20, 65, 60, 20])
plot(0:pi/1e3:pi, evmiuw1(1:1001))
plot(0:pi/1e3:pi, evmiuw(1:1001))
set(but1, 'position', [10, 65, 60, 20])
set(but1, 'position', [0, 65, 60, 20])
set(but, 'position', [0, 0, 60, 20])
set(but1, 'position', [0, 70, 60, 20])
set(but1, 'position', [70, 0, 60, 20])
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
tic; drawnow; toc
stop
edit dctIfrom_ig
edit dctI
edit dctIintgrid
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
figure
plot(0:pi/1e3:pi, evmiuw(1:1001))
plot(0:pi/1e3:pi, fieldw(1:1001))
plot(0:0.2:1e3, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
max(imag(allevmiut))
max(imag(allevmiut)./real(allevmiut))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
J1
figure
plot(0:pi/1e3:pi, evmiuw(1:1001))
max(imag(evmiuw))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0, @(x) 0, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 10);
figure
plot(0:pi/1e3:pi, evmiuw(1:1001))
hold on
plot(0:pi/1e3:pi, evmiuw(1:1001), 'r')
load gong.mat;
gong = audioplayer(y, Fs);
play(gong);
whos
Fs
clear y Fs
figure
xdata = 0:0.01:1;
plot(xdata, 0.5*(tanh(50*(xdata-0.9)) - tanh(5))
plot(xdata, 0.5*(tanh(50*(xdata-0.9)) - tanh(5)))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.8)) - tanh(10)), @(x) 25*sech(50*(x-0.8))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
weight
nprop
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.5)) - tanh(25)), @(x) 25*sech(50*(x-0.5))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
whos
clear all
load optV40
whos
Vf = @(x)1-1./sqrt(x.^2+1)
xdomain = [-120 120];
xdlength = 240;
Nx = 128*3
dx=xdlength/Nx
[Vabs, xabs, x, p, K, Nx] = get_prop_vars(Vf, xdomain, dx, Vopt649x);
[fi0, E0, x, E, P, H] = gsV(Vabs, xdomain, Nx);
E(1:10)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
norm(psi(:, end))
figure
plot(0:0.2:1e3, fieldt)
viewVP(psi, Vf, fieldt, x, 0.01)
figure
plot(0:pi/1e3:pi, evmiuw(1:1001), 'r')
J1
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.5*(1-tanh(100*(w-0.07))), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
J1
norm(psi(:, end))
J11
figure
plot(0:0.2:1e3, fieldt)
plot(0:0.2:1e3, fieldt1)
plot(0:pi/1e3:pi, evmiuw1(1:1001), 'r')
viewVP(psi1, Vf, fieldt1, x, 0.01)
whos
save xdomain120
%-- 20/02/2014 14:01 --%
load coulomb_optV40
w=0:pi/1e3:pi/0.2;
fieldwr = (rand(1, 5001)-0.5).*0.5.*(1-tanh(100*(w-0.07)));
plot(w(1:101), fieldwg(1:101))
plot(w(1:101), fieldwr(1:101))
fieldtg = dctI(fieldwr)/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtr = dctI(fieldwr)/dctfactor;
t=0:0.2:1e3;
plot(t, filedtg)
plot(t, filedtr)
plot(t, fieldtr)
[fieldt3, fieldw3, psi3, evmiut3, evmiuw3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, weight3] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwg, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt3, fieldw3, psi3, evmiut3, evmiuw3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, weight3] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwr, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-5, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
[fieldt3, fieldw3, psi3, evmiut3, evmiuw3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, weight3] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwr, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
weight
sin(eps)/eps
sin(5*eps)/eps
norm(allpsi(:, end))
weight
figure
norm(psi(:, end))
norm(psi3(:, end))
J1
J13
figure
plot(w, fieldw)
plot(w, fieldw3)
plot(w(1:101), fieldw3(1:101))
figure
plot(t, fieldt)
plot(t, fieldt3)
figure
plot(w(1:1001), evmiuw(1:1001))
plot(w(1:1001), evmiuw3(1:1001))
figure
plot(t, evmiut3)
viewVPmiux(psi3, Vf, fieldt3, xabs, 0.01)
viewVPmiux(psi3, Vf, xabs, fieldt3, x, 0.01)
viewVPmiux(psi3, Vabs, xabs, fieldt3, x, 0.01)
fieldwr1 = (rand(1, 5001)-0.5).*0.5.*(1-tanh(100*(w-0.07)));
fieldtr1 = dctI(fieldwr1)/dctfactor;
figure
plot(t, fieldtr)
plot(t, fieldtr1)
[fieldt4, fieldw4, psi4, evmiut4, evmiuw4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, weight4] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwr1, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(psi4(:, end))
figure
plot(w(1:101), fieldw4(1:101))
figure
J14
figure
plot(t, fieldt4)
figure
plot(w(1:1001), evmiuw4(1:1001))
[fieldt4, fieldw4, psi4, evmiut4, evmiuw4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, weight4] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, zeros(1, 5001), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt4, fieldw4, psi4, evmiut4, evmiuw4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, weight4] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwr*1e-7, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt4, fieldw4, psi4, evmiut4, evmiuw4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, weight4] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwr*1e-7, @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(psi4(:, end))
J14
%-- 23/02/2014 14:28 --%
load coulomb_optV40
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sum(J1fun)
sum(J2fun)
penalnormf(normpsiT)
normpsiT
weight
conv(70:80)
conv(70:82)
sum(J1fun)
sum(J2fun)
normpsiT
penalnormf(normpsiT)
norm(psi(:, end))
J1
conv(end)
figure
plot(w(1:101), fieldw(1:101))
w=0:pi/1e3:pi/0.2;
t=0:0.2:1e3;
plot(t, fieldt)
figure
plot(w, fieldw)
plot(w(1:101), fieldw(1:101))
plot(w(1:1001), evmiuw(1:1001))
fieldiw = instw(fieldt, t);
figure
plot(t, fieldtiw)
plot(t, fieldiw)
doc hilbert
fieldiw = instw(hilbert(fieldt), t);
plot(t, fieldiw)
fieldiw = instw(fieldt, 0.2);
plot(t, fieldiw)
fieldiw = instw(hilbert(fieldt), t);
plot(t, fieldiw)
fieldiw = instw(hilbert(fieldt), dt);
fieldiw = instw(hilbert(fieldt), 0.2);
plot(t, fieldiw)
fieldiw1 = instw(fieldt(1:5:end), 1);
figure
plot(t(1:5:end), fieldiw1)
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(evmiut(1:5:end), 1, 250);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(fieldt(1:5:end), 1, 250);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(fieldt(1:10:end), 2, 125);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(fieldt(1:20:end), 4, 62);
[wreal, d_result, w_result, g_result, meand, w_maxd, maxd] = KBDMrinp(fieldt(1:40:end), 8, 31);
figure
plot(w(1:101), fieldw(1:101))
figure
plot(wreal, meand)
plot(wreal(meand<1e2), meand(meand<1e2))
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e1*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-1, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sum(J1fun)
sum(J2fun)
penalnormf(normpsiT)
max(abs(allfield))
max(abs(fieldg))
plot(w, fguess(w))
plot(0:dw:pi/T, fguess(w))
plot(0:dw:pi/T, fguess(0:dw:pi/T))
norm(psi(:, end))
J11
J1
figure
plot(t, fieldt1)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, weight, 1e3, 0.2, 7, 7, 1e-4, 1e3
fieldtg = dctI(0.6*0.5*(1-tanh(100*(w-0.07))))/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(0.6*0.5*(1-tanh(100*(w-0.07))))/dctfactor;
figure
plot(t, fieldtg
plot(t, fieldtg)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, weight, 1e3, 0.2, 7, 7, 1e-4, 1e3);
norm(allpsi(:, end))
relE
norm(psi(:, end))
J1
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, weight, 1e3, 0.2, 7, 7, 1e-4, 2e3
weight
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, weight, 1e3, 0.2, 7, 7, 1e-4, 2e3);
relE
weight
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, weight/2, 1e3, 0.2, 7, 7, 1e-4, 2e3);
J1
norm(psi(:, end))
clear psi1
save tol4
%-- 24/02/2014 13:48 --%
load tol4
whos
figure
plot(w(1:101), fieldw(1:101))
plot(t, fieldt)
figure
plot(t, evmiut)
plot(w(1:1001), evmiuw(1:1001))
2*pi/85.2
figure
plot(w(1:1001), evmiuw(1:1001))
plot(w(1:101), fieldw(1:101))
2*pi/84.6
2*pi/(168.7-84.6)
2*pi/(251.6-168.7)
2*pi/(716.4-512.4)
2*pi/(944.6-716.4)
E(2:10)-E(1:9)
whos
psiE = P\psi;
figure
plot(t, conj(psiE).*psiE)
viewVPmiux(psi, Vabs, xabs, fieldt, x, 0.01)
plot(t, conj(psiE(2:end)).*psiE(2:end))
plot(t, conj(psiE(:,2:end)).*psiE(:,2:end))
plot(t, conj(psiE(2:end, :)).*psiE(2:end, :))
figure
plot(t, evmiut)
fieldwa=dctI(fieldt(1:(251.8/0.2 +1)));
figure
plot(w(1:101), fieldwa(1:101))
251.8/0.2
plot(0:pi/251.8:100*pi/251.8, fieldwa(1:101))
fieldwb=dctI(fieldt((293.6/0.2 + 1):end));
figure
1e3-293.6
plot(0:pi/706.4:100*pi/706.4, fieldwb(1:101))
fieldiw = instw(fieldt(1:5:end), 1);
figure
plot(t(1:5:end), fieldiw)
isreal(fieldt)
phi=angle(hilbert(fieldt));
phi = phi_cont(phi);
fieldiw2 = (phi(2:end)-phi(1:end-1))/0.2;
figure
plot(t(2:end), fieldiw2)
fieldiw2 = (phi(2:5:end)-phi(1:5:end-1));
plot(t(2:5:end), fieldiw2)
fieldiw2 = (phi(2:end)-phi(1:end-1))/0.2;
plot(t(2:end), fieldiw2)
fieldiw = instw(fieldt(1:5:end), 1);
figure
plot(0:1e3, phi)
figure
plot(0:1e3, w)
phicont = phi_cont(phi);
figure
plot(0:1e3, phicont)
w(1) = (phi(2) - phi(1))./dx(1);
w(2:(N - 1)) = (phi(3:N) - phi(1:(N - 2)))./(dx(1:(N - 2)) + dx(2:(N - 1)));
w(N) = (phi(N) - phi(N - 1))./dx(N - 1);
figure
plot(0:1e3, w)
wc(1) = (phicont(2) - phicont(1))./dx(1);
wc(2:(N - 1)) = (phicont(3:N) - phicont(1:(N - 2)))./(dx(1:(N - 2)) + dx(2:(N - 1)));
wc(N) = (phicont(N) - phicont(N - 1))./dx(N - 1);
plot(0:1e3, wc)
fieldiw = instw(fieldt, 0.2);
plot(t, fieldiw)
fieldiw = instw(fieldt(1:5:end), 1);
plot(t(1:5:end), fieldiw)
fieldiw = instw(fieldt, 0.2);
figure
plot(0:0.2:1e3, w)
mod(-eps, 3)
mod(eps, 3)
fieldiw = instw(fieldt, 0.2);
plot(0:0.2:1e3, fieldiw)
save tol4
%-- 26/02/2014 12:47 --%
load tol4
figure
plot(t, fieldiw)
W = psi2wigner([fieldt(end:-1:2), fieldt(1:end-1)]);
%-- 26/02/2014 13:55 --%
load tol4
clear psiE
save tol4
W = psi2wigner([fieldt(end:-10:2), fieldt(1:10:end-1)]);
W = psi2wigner([fieldt(end:-10:2), fieldt(1:10:end-1)].');
mesh(W)
xint = xdomain(1):dx/2:(xdomain(2) - dx/2);
tW = -1e3:(1e3 - 1);
wW = -pi/2:pi/2e3:(pi/2 - pi/2e3);
mesh( tW, wW, W)
mesh(tW(1001:2000), wW, W(1001:2000, :))
mesh(tW(1001:2000), wW, W(:, 1001:2000))
figure
plot(wW, W(:,1100))
plot(wW(1001:2000), W(1001:2000,1100))
plot(wW(1001:1300), W(1001:1300,1100))
plot(wW(1001:1100), W(1001:1100,1100))
plot(wW(1001:1100), W(1001:1100,1001))
plot(wW(1001:1100), W(1001:1100,1101))
plot(wW(1001:1100), W(1001:1100,1201))
plot(wW(1001:1100), W(1001:1100,1301))
figure
plot(t, fieldt)
plot(wW(1001:1100), W(1001:1100,1401))
plot(wW(1001:1100), W(1001:1100,1501))
plot(wW(1001:1100), W(1001:1100,1601))
plot(wW(1001:1100), W(1001:1100,1701))
plot(wW(1001:1100), W(1001:1100,1710))
plot(wW(1001:1100), W(1001:1100,1720))
plot(wW(1001:1100), W(1001:1100,1820))
plot(wW(1001:1100), W(1001:1100,1920))
plot(wW(1001:1100), W(1001:1100,1300))
plot(wW(1001:1100), W(1001:1100,1200))
plot(wW(1001:1100), W(1001:1100,1100))
plot(wW(1001:1100), W(1001:1100,1500))
plot(wW(1001:1100), W(1001:1100,1600))
figure
plot(t, fieldiw)
plot(wW(1001:1100), W(1001:1122,1600))
plot(wW(1001:1100), W(1001:1100,1122))
plot(wW(1001:1100), W(1001:1100,1063))
figure
plot(0:pi/251.8:100*pi/251.8, fieldwa(1:101))
plot(wW(1001:1100), W(1001:1100,1590))
figure
plot(0:pi/706.4:100*pi/706.4, fieldwb(1:101))
plot(wW(1001:1100), W(1001:1100,1800))
figure
plot(t, fieldt)
plot(wW(1001:1100), W(1001:1100,1700))
mesh(tW(1001:2000), wW(1001:2000), W(1001:2000, 1001:2000))
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 20*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(allpsi(:, end))
penalnormf(normpsiT)
sum(J1fun)
conv(2)
norm(allpsi(:, end))
penalnormf(normpsiT)
penalnormf(9.94)
penalnormf(9.9)
penalnormf(1)
penalnormf(0.994)
penalnormf(0.99)
penalnormf(0.98)
sum(J1fun)
clear psi1
[fieldt2, fieldw2, psi2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, weight2] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.65, 0.67), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
chebw = chebweights_gen(5, [0 1], [0:0.1:1])
(-2*result_points + domain(1) + domain(2))/lengthD
chebw = chebweights_gen(5, [0 1], [0:0.1:1])
theta = 0:pi/4:pi
y=cos(theta)
z = y/2 + 1
z = (y+1)/2
z = (y(end:-1:1)+1)/2
fz = z.^2
fz*chebw
z = (y+1)/2
fz = z.^2
fz*chebw
(0:0.1:1).^2/3
(0:0.1:1).^3/3
z = (y(end:-1:1)+1)/2
fz = z.^2
chebw = chebweights_gen(5, [0 1], [0:0.1:1])
z = (y+1)/2
z = (y(end:-1:1)+1)/2
fz = z.^2
fz*chebw
(0:0.1:1).^3/3
(z.^4)*chebw
z.^5/5
(0:0.1:1).^5/5
(z.^5)*chebw
(0:0.1:1).^6/6
z = (y(end:-1:1)+1)
chebw = chebweights_gen(6, [0 2], [0:0.2:2])
z = (y(end:-1:1)+1)
(z.^5)*chebw
theta = 0:pi/5:pi
y=cos(theta)
z = (y(end:-1:1)+1)
(z.^5)*chebw
z.^6/6
(0:0.1:1).^6/6
(0:0.2:2).^6/6
(z.^5)*chebw-(0:0.2:2).^6/6
sum(chebw)
%-- 03/03/2014 19:59 --%
%-- 12/03/2014 15:33 --%
whos
relE
figure
plot(0:nitertot, convtot)
figure
w=0:pi/2e3:pi/0.2;
t=0:0.2:2e3;
plot(w, fieldw)
plot(w(1:101), fieldw(1:101))
figure
plot(t, fieldt)
psiE=P\psi;
figure
plot(t, conj(psiE).*psiE)
plot(t, conj(psiE(1:5, :)).*psiE(1:5, :))
plot(t, conj(psiE(2:7, :)).*psiE(2:7, :))
clear all
whos
relE
figure
plot(t, fieldt)
figure
plot(w(1:1001), fieldw(1:1001))
hold on
plot(w(1:1001), evmiuw(1:1001))
doc plotyy
figure
DocPropFig
figure
DocPropFig
set(AX(2), 'Ylim', [-20 35])
set(AX(2), 'Ylim', [-30 35])
set(AX(1), 'Ylim', [-30 35]/30)
set(AX(2), 'Ylim', [-20 35])
set(AX(1), 'Ylim', [-20 35]/20)
set(AX(1), 'Ylim', [-20 35]/15)
set(AX(2), 'Ylim', [-20 35]/35*1.8)
set(AX(2), 'Ylim', [-20 35])
set(AX(1), 'Ylim', [-20 35]/35*1.8)
set(AX(2), 'Ytick', [-20, -10, 0, 10, 20, 30])
set(AX(1), 'Ytick', [-1, -0.5, 0, 0.5, 1, 1.5])
set(AX(1), 'Xlim', [0 0.6])
set(AX(2), 'Xlim', [0 0.6])
set(get(AX(2),'YLabel'), 'Position', [0, 0.28, 0])
figure
xlabel('$t$ (a.u.)', 'interpreter', 'latex')
%-- 13/03/2014 10:41 --%
load('C:\Users\Ido\Dropbox\MATLAB\OCT\HHGcontrol\tol4.mat')
ylabel('$\epsilon(t)$ (a.u.)', 'interpreter', 'latex')
norm(psi(:, end))
%-- 19/03/2014 11:30 --%
load('C:\Users\Ido\Dropbox\MATLAB\OCT\HHGcontrol\tol4.mat')
whos
size(fieldw)
fieldrand = rand(1, 5001);
max(abs(fieldw))
fieldrand = rand(1, 5001).*0.6*0.5*(1-tanh(100*(w-0.07)));
fieldrand = rand(1, 5001).*0.6*0.5*(1-tanh(100*(w.'-0.07)));
figure
plot(w(1:101), fieldrand(1:101))
size(w)
size(fieldrand)
fieldrand = rand(1, 5001)*0.6*0.5.*(1-tanh(100*(w-0.07)));
plot(w(1:101), fieldrand(1:101))
fieldrand = (rand(1, 5001)-0.5)*0.6*0.5.*(1-tanh(100*(w-0.07)));
fieldrand = (rand(1, 5001)-0.5)*0.6.*(1-tanh(100*(w-0.07)));
plot(w(1:101), fieldrand(1:101))
fieldrand = (rand(1, 5001)-0.5).*(1-tanh(100*(w-0.07)));
plot(w(1:101), fieldrand(1:101))
fieldrandt = dctI(fieldrand)/dctfactor;
figure
plot(t, fieldrandt)
plot(t, fieldt)
hold on
plot(t, fieldrandt, 'r')
fieldwp = fieldw + 1e-2*fieldrand;
figure
plot(w(1:101), fieldw(1:101))
hold on
plot(w(1:101), fieldrand(1:101), 'r')
plot(w(1:101), fieldwp(1:101), 'r')
fieldtp = dctI(fieldwp)/dctfactor;
plot(t, fieldtp, 'g')
[fieldtp0, fieldwp0, psip, evmiutp, evmiuwp, mnitercp, Jp, J1p, J2p, Jorthp, Jpnormp] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwp, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
max(abs(fieldtp0-fieldtp))
max(abs(fieldtwp0-fieldwp))
max(abs(fieldwp0-fieldwp))
J
conv(end)
J
Jp
J1
J1p
norm(ps
norm(psi(:, end))
norm(psip(:, end))
Jpnorm
Jpnormp
J2
J2p
whos
conv1(end)
Jp
Jorthp
J1p+J2p+Jpnormp
J1p
J1
[fieldtp0, fieldwp0, psip, evmiutp, evmiuwp, mnitercp, Jp, J1p, J2p, Jorthp, Jpnormp] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwp, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
Jp
conv(end)
fieldwp = fieldw + 5e-2*fieldrand;
figure
fieldtp = dctI(fieldwp)/dctfactor;
plot(t, fieldtp, 'g')
[fieldtp0, fieldwp0, psip, evmiutp, evmiuwp, mnitercp, Jp, J1p, J2p, Jorthp, Jpnormp] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwp, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
Jp
figure
plot(w, evmiuw)
plot(w(1:1001), evmiuw(1:1001))
hold on
plot(w(1:1001), evmiuwp(1:1001), 'r')
J1p
J1
fieldwp2 = fieldw + 1e-1*fieldrand;
fieldtp2 = dctI(fieldwp2)/dctfactor;
plot(t, fieldtp, 'm')
[fieldtp20, fieldwp20, psip2, evmiutp2, evmiuwp2, mnitercp2, Jp2, J1p2, J2p2, Jorthp2, Jpnormp2] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwp, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldtp20, fieldwp20, psip2, evmiutp2, evmiuwp2, mnitercp2, Jp2, J1p2, J2p2, Jorthp2, Jpnormp2] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwp2, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
Jp2
J1p2
conv(end)
plot(w(1:1001), evmiuwp2(1:1001), 'g')
E(2:10)-E(1)
whos
psiE=P\psi;
figure
plot(t, conj(psiE(2:7, :)).*psiE(2:7, :))
E(3:10)-E(2)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.44, 0.46), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
clear all
load coulomb_optV40
whos
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.44, 0.46), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.47, 0.49), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
figure
w=0:pi/1e3:pi/0.2;
t=0:0.2:1e3;
plot(w(1;1001), evmiuw(1:1001))
plot(w(1:1001), evmiuw(1:1001))
figure
plot(w(1:101), fieldw(1:101))
figure
plot(t, fieldt)
figure
plot(t, evmiut)
psiE=P\psi;
figure
plot(t, conj(psiE).*psiE)
viewVPmiux(psi, Vabs, xabs, fieldt, x, 0.01)
norm(psi(:, end))
save field48 fieldt fieldw evmiut relE conv niter mallniterc J1 maxgrad weight
fieldiw = instw(fieldt, 0.2);
figure
plot(t, fieldiw)
W = psi2wigner([fieldt(end:-10:2), fieldt(1:10:end-1)].');
mesh(W)
tW = -1e3:(1e3 - 1);
wW = -pi/2:pi/2e3:(pi/2 - pi/2e3);
mesh( tW, wW, W)
mesh(tW(1001:end), wW(1001:end), W)
mesh(tW(1001:end), wW(1001:end), W(1001:end, 1001:end))
mesh(tW(1001:end), wW(1001:1101), W(1001:end, 1101:end))
mesh(tW(1001:end), wW(1001:1101), W(1001:end, 1001:1101))
size(W)
figure
plot(wW(1001:1101
plot(wW(1001:1101), W(1001, 1001:1101))
plot(wW(1001:1501), W(1001, 1001:1501))
plot(wW(1001:end), W(1001, 1001:end))
plot(wW(1001:end), W(1101, 1001:end))
plot(wW(1001:end), W(1201, 1001:end))
plot(wW(1001:end), W(1301, 1001:end))
clear W
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.47, 0.49), 0, 1e-2, 2e3, 0.2, 7, 7, 1e-3, 1e3);
figure
w=0:pi/2e3:pi/0.2;
t=0:0.2:2e3;
plot(t, fieldt)
norm(psi(:, end))
figure
plot(w, fieldw)
plot(w(1:101), fieldw(1:101))
figure
plot(w(1:1001), evmiuw(1:1001))
save field48_2000 fieldt fieldw evmiut relE conv niter mallniterc J1 maxgrad weight
save field48_2000 fieldt fieldw evmiut relE conv niter mallniterc J1 maxgrad weight w t
%-- 24/03/2014 09:22 --%
load coulomb_optV40
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
norm(psi(:, end))
viewVPmiux(psi, Vabs, xabs, fieldt, x, 0.01)
figure
whos
w=0:pi/1e3:pi/0.2;
t=0:0.2:1e3;
plot(w, fieldw)
plot(w(1:101), fieldw(1:101))
plot(w(1:1001), evmiuw(1:1001))
figure
plot(t, fieldt)
psiE=P\psi;
figure
plot(t, conj(psiE).*psiE)
E(2:10)-E(1)
E(3:10)-E(2)
plot(w, evmiuw)
pl0t(w(1:1001
plot(w(1:1001), evmiuw(1:1001))
E(4:10)-E(3)
save field30 fieldt fieldw evmiut relE conv niter mallniterc J1 maxgrad weight w t
J1
%-- 25/03/2014 09:20 --%
load('field30.mat')
load coulomb_optV40
[~, ~, psi, ~, ~, ~, J, ~, ~, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
J
conv(end)
Jpnorm
figure
plot(t, fieldt)
2*pi/100
2*pi/(462.8-357.4)
2*pi/(562.4-462.8)
figure
plot(w(1:101), fieldw(1:101))
dctfactor = 1e3/(sqrt(5e3*pi))
0.06/pi*1e3
19*pi/1e3
20*pi/1e3
19*pi/1e3*5
fieldwf = zeros(1, 1001);
fieldwf(20)=40;
fieldft=dctI(fieldwf)/dctfactor;
figure
fieldwf = zeros(1, 5001);
fieldwf(20)=40;
fieldft=dctI(fieldwf)/dctfactor;
plot(t, fieldtf)
fieldtf=dctI(fieldwf)/dctfactor;
clear fieldft
plot(t, fieldtf)
fieldwf(20)=20;
fieldtf=dctI(fieldwf)/dctfactor;
plot(t, fieldtf)
[fieldtf, fieldwf, psif, evmiutf, evmiuwf, mnitercf, Jf, J1f, J2f, Jorthf, Jpnormf] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldwf, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
J
Jf
Jpnormf
viewVPmiux(psi, Vabs, xabs, fieldtf, x, 0.01)
J1f
norm(psif(:, end))
J2f
figure
plot(w, evmiuwf)
plot(w(1:101), evmiuwf(1:101))
plot(w(1:1001), evmiuwf(1:1001))
figure
plot(w(1:1001), evmiuw(1:1001))
whos
evmiuw = dctI(evmiut)*dctfacotr;
evmiuw = dctI(evmiut)*dctfactor;
[~, ~, psi, ~, evmiuw, ~, J, ~, ~, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
save field30 evmiuw -append
load('field30.mat')
figure
plot(w, evmiuw)
plot(w(1:1001), evmiuw(1:1001))
psiEf=P\psif;
figure
plot(t, conj(psiEf).*psiEf)
viewVPmiux(log(norm(psi)), Vabs, xabs, fieldtf, x, 0.01)
viewVPmiux(log(abs(psi)), Vabs, xabs, fieldtf, x, 0.01)
viewVPmiux(psif, Vabs, xabs, fieldtf, x, 0.01)
norm(psif(:, end))
viewVPmiux(abs(psif), Vabs, xabs, fieldtf, x, 0.01)
viewVPmiux(log10(abs(psif)), Vabs, xabs, fieldtf, x, 0.01)
figure
plot(x, U(:, 1))
plot(x, P(:, 1))
viewVlogPmiux(psif, Vabs, xabs, fieldtf, x, 0.01)
whos
[Vabs160, xabs160, x160, p160, K160, Nx160, V0160] = get_prop_vars(Vf, [-160 160], dx, Vopt649x);
figure
plot(x160, real(Vabs160))
Nx160
xdomain
xdomain160=[-160 160];
save coulomb_optV120 Vabs160 xabs160 x160 p160 K160 Nx160 V0160 Vopt649x
[fieldt120, fieldw120, psi120, evmiut120, evmiuw120, mniterc120, J120, J1120, J2120, Jorth120, Jpnorm120] = guessresultspnorm(fi0160, Vabs160, 1, xdomain160, xabs160, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
[fi0160, E0160, x160, E160, P160, H160] = gsV(Vf, xdomain160, Nx160);
save coulomb_optV120 Vabs160 xabs160 x160 p160 K160 Nx160 V0160 Vopt649x fi0160 E0160 x160 E160 P160 H160
[fieldt120, fieldw120, psi120, evmiut120, evmiuw120, mniterc120, J120, J1120, J2120, Jorth120, Jpnorm120] = guessresultspnorm(fi0160, Vabs160, 1, xdomain160, xabs160, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
J120
J
norm(psi120(:, end))
norm(psi(:, end))
Jpnorm120
Jpnorm
J1120
J2120
J2
[~, ~, ~, ~, ~, ~, ~, ~, J2, ~, ~] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
J2
J2120
figure
plot(w(1:1001), evmiuw120(1:1001))
figure
plot(t, evmiut)
hold on
plot(t, evmiut120, 'r')
viewVPmiux(psi120, Vabs120, xabs120, fieldt120, x120, 0.01)
viewVPmiux(psi120, Vabs160, xabs160, fieldt120, x120, 0.01)
viewVPmiux(psi120, Vabs160, xabs160, fieldt120, x160, 0.01)
figure
plot(x, xabs)
112.8/0.2
evmiu113 = evmiu(psi(:,565), xabs, 1);
figure
max(abs(evmiu113-evmiut(565)))
evmiu120_113 = evmiu(psi120(:,565), xabs120, 1);
evmiu120_113 = evmiu(psi120(:,565), xabs160, 1);
max(abs(evmiu120_113-evmiut120(565)))
figure
plot(x, conj(psi(:,565).*psi(:,565)))
plot(x, conj(psi(:,565)).*psi(:,565))
hold on
plot(x160, conj(psi120(:,565)).*psi120(:,565), 'r')
figure
plot(x, conj(psi120(65:(512-64),565)).*psi120(:,565), 'r')
512-64
plot(x, conj(psi120(65:448,565)).*psi120(65:448,565), 'r')
plot(x, conj(psi120(129:384,565)).*psi120(129:384,565), 'r')
plot(x, conj(psi120(129:384,565)).*psi120(129:384,565)-conj(psi(:,565)).*psi(:,565), 'r')
evmiu120_80 = evmiu(psi120(129:384,:), xabs, 1);
figure
plot(t, evmiu120_80)
plot(t, evmiu120_80, 'g')
hold on
plot(w(1:1001), evmiuw(1:1001), 'r')
figure
plot(w(1:1001), evmiuw120(1:1001)-evmiuw(1:1001), 'r')
plot(w(1:1001), (evmiuw120(1:1001)-evmiuw(1:1001))./evmiuw120(1:1001))
evmiu120_80_2 = evmiu(psi120(129:384,:), xabs160(129:384), 1);
figure
plot(t, evmiu120_80_2, 'm')
evmiu120_120 = evmiu(psi120(65:448,:), xabs160(65:448), 1);
figure
plot(t, evmiu120_120, 'k')
[Vabs240, xabs240, x240, p240, K240, Nx240, V0240] = get_prop_vars(Vf, [-200 200], dx, Vopt649x);
[fi0240, E0240, x240, E240, P240, H240] = gsV(Vf, xdomain240, Nx240);
xdomain160=[-240 240];
[fi0240, E0240, x240, E240, P240, H240] = gsV(Vf, xdomain240, Nx240);
xdomain160=[-160 160];
xdomain240=[-240 240];
[fi0240, E0240, x240, E240, P240, H240] = gsV(Vf, xdomain240, Nx240);
[fieldt240, fieldw240, psi240, evmiut240, evmiuw240, mniterc240, J240, J1240, J2240, Jorth240, Jpnorm240] = guessresultspnorm(fi0240, Vabs240, 1, xdomain240, xabs240, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
J240
J120
J1240
Jpnorm120
Jpnorm240
norm(psi(:, end))
norm(psi120(:, end))
norm(psi240(:, end))
plot(t, evmiut240, 'c')
viewVPmiux(psi240, Vabs240, xabs240, fieldt240, x240, 0.01)
plot(w(1:1001), evmiuw240(1:1001), 'g')
E240(1:10)
[E(1:10) E240(1:10)]
[E(1:10)-E(1) E240(1:10)-E240(1)]
mniterc240
mniterc
mnitercf
mniterc120
max(abs(fieldt240-fieldt))
[Vabs240, xabs240, x240, p240, K240, Nx240, V0240] = get_prop_vars(Vf, [-240 240], dx, Vopt649x);
xdomain240=[-240 240];
[fi0240, E0240, x240, E240, P240, H240] = gsV(Vf, xdomain240, Nx240);
[fieldt240, fieldw240, psi240, evmiut240, evmiuw240, mniterc240, J240, J1240, J2240, Jorth240, Jpnorm240] = guessresultspnorm(fi0240, Vabs240, 1, xdomain240, xabs240, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.29, 0.31), 0, 1e3, 0.2, 7, 7, 1e-3);
J240
J120
J1240
J1120
Jpnorm240
norm(psi240(:, end))
norm(psi(:, end))
plot(t, evmiut240, 'c')
viewVPmiux(psi240, Vabs240, xabs240, fieldt240, x240, 0.01)
plot(w(1:1001), evmiuw240(1:1001), 'g')
save coulomb_optV240 Vabs240 xabs240 x240 p240 K240 Nx240 V0240 Vopt649x fi0240 E0240 x240 E240
whos
clear psi psi120 psi240 psiEf psif H H160 H240 P P160 P240
whos
save
%-- 08/05/2014 13:38 --%
[fieldtp, fieldwp, psip, evpt, evpw, evmiutp, evmiuwp, relEp, convp, niterp, mallnitercp, J1p, maxgradp, weightp] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.47, 0.49)./(w.^2+eps), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
load coulomb_optV40
[fieldtp, fieldwp, psip, evpt, evpw, evmiutp, evmiuwp, relEp, convp, niterp, mallnitercp, J1p, maxgradp, weightp] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.47, 0.49)./(w.^2+eps), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtp, fieldwp, psip, evpt, evpw, evmiutp, evmiuwp, relEp, convp, niterp, mallnitercp, J1p, maxgradp, weightp] = OCfpnorm_evp(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.47, 0.49)./(w.^2+eps), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(w, evpw)
w=0:pi/1e3:pi/0.2;
b=89
plot(0:pi/1e3:pi, evp(1:1001))
plot(0:pi/1e3:pi, evpw(1:1001))
figure
plot(0:0.2:1e3, evpt)
plot(0:0.2:1e3, evpallt(1:6:end))
plot(0:0.2:1e3, allevpt(1:6:end))
figure
plot(0:pi/1e3:pi, evpw(1:1001))
figure
plot(0:0.2:1e3, allevpt(1:6:end))
relE
figure
plot(0:0.2:1e3, fieldt)
plot(0:pi/1e3:pi/10, fieldw(1:101))
fieldt = dctI(fieldw)/dctfactor;
figure
plot(0:0.2:1e3, fieldt)
figure
plot(0:pi/1e3:pi, evpw(1:1001))
figure
plot(0:pi/1e3:pi, evpw(1:1001))
figure
plot(0:pi/1e3:pi, evmiuwp(1:1001))
figure
plot(0:pi/1e3:pi, evpw(1:1001)./(0:pi/1e3:pi))
norm(psi(:, end))
norm(psip(:, end))
figure
plot(0:0.2:1e3, fieldtp)
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(w(1:101), fieldt(1:101))
plot(w(1:101), fieldtp(1:101))
plot(w(1:101), fieldwp(1:101))
figure
plot(0:niter, conv)
figure
plot(w(1:101), fieldw(1:101))
figure
plot(t, fieldt)
figure
plot(w(1:1001), evmiuw)
whos
evmiuw = dctI(evmiut)*dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
evmiuw = dctI(evmiut)*dctfactor;
figure
plot(w(1:1001), evmiuw)
plot(w(1:1001), evmiuw(1:1001))
psiEp=P\psip;
figure
plot(t, conj(psiEp).*psiEp)
E(2:10)-E(1:9)
E(3:10)-E(1:8)
E(4:10)-E(1:7)
E(4:10)-E(2:8)
E(5:10)-E(1:6)
E(6:10)-E(1:5)
E(7:10)-E(1:4)
E(8:10)-E(1:3)
E(9:10)-E(1:2)
save field48p fieldtp fieldwp evmiutp evmiuwp relEp convp niterp mallnitercp J1p maxgradp weightp evpt evwp
save field48p fieldtp fieldwp evmiutp evmiuwp relEp convp niterp mallnitercp J1p maxgradp weightp evpt evpw
fieldiw = instw(fieldt);
fieldiw = instw(fieldt, 0.2);
figure
plot(t, fieldiw)
fieldiwp = instw(fieldtp, 0.2);
figure
plot(t, fieldiwp)
%-- 11/05/2014 18:06 --%
%-- 13/05/2014 10:51 --%
load coulomb_optV40
%-- 13/05/2014 11:04 --%
pwd
userpath('C:\Users\Ido\Dropbox\MATLAB')
load coulomb_optV40
%-- 13/05/2014 11:06 --%
pwd
path
%-- 13/05/2014 11:13 --%
load coulomb_optV40
67.9
phi = angle(hilbert(fi0));
figure
plot(x, phi)
phi(1) = phi(2);
plot(x, phi)
iwfi0 = Dcosines(fi0, 160);
figure
plot(x, iwfi0)
iwfi0 = Dcosines(phi, 160);
plot(x, iwfi0)
[xdist, xint] = xdistW(fi0, xdomain, @(p) p);
figure
plot(xint, xdist)
[xdist, xint] = xdistW(real(fi0), xdomain, @(p) p);
plot(xint, xdist)
[xdist, xint] = xdistW(real(fi0), xdomain, @(p) abs(p));
figure
plot(xint, xdist)
figure
plot(x, iwfi0)
sum(xdist)
mp = evp(fi0, p);
p(1:10)
whos
figure
fi0p = fft(fi0);
sum(fi0p'.*fi0p)
fi0p'*fi0p
ans*1/sqrt(2*pi)
fi0p'*fi0p/256
size(fi0p)
fi0p = fft(fi0)/sqrt(256);
fi0p'*fi0p
figure
plot(p, fi0p)
plot(p, conj(fi0p).*fi0p)
fi0p'*(p.*fi0p)
fi0p'*(abs(p).*fi0p)
mean(fi0iw)
mean(iwfi0)
viewWigner(fi0, xdomain, 1, [0 0.1])
Wrange = getWrange(fi0)
viewWigner(fi0, xdomain, 1, Wrange)
pi/0.6
pi/sqrt(2)
figure
plot(x, fi0)
iwV = Dcosines(V0, 160);
whos
iwV = Dcosines(Vf(x), 160);
figure
plot(x, iwV)
phiV = Dcosines(Vf(x), 160);
phiV = phi_cont(angle(hilbert(Vf(x))));
figure
plot(x, phiV)
iwV = Dcosines(phiV, 160);
figure
plot(x, phiV)
plot(x, iwV)
pi/4.2
iwV = Dcosines([phiV; phiV(end)], 160);
plot(x, iwV)
plot(x, iwV(1:end-1))
whos
save iw_mapping xint xdist phiV phi iwfi0 iwV fi0p Wrange
load tol4
iw1000 = instwcos(psi(:, end));
iw1000 = instwcos(psi(:, end), 160);
figure
plot(x, psi(:,end))
plot(x, imag(psi(:,end)))
hold on
plot(x, psi(:,end), 'r')
iwr1000 = instwcos(real(psi(:, end)), 160);
figure
plot(x, irw1000)
plot(x, iwr1000)
iwi1000 = instwcos(imag(psi(:, end)), 160);
figure
plot(x, iwi1000)
iw1000 = instwcos(psi(:, end));
iw1000 = instwcos(psi(:, end), 160);
figure
plot(x, iw1000)
%-- 14/05/2014 12:24 --%
load tol4
evpt = evp(psi, p);
dctfactor = 1e3/(sqrt(5e3*pi))
evpw = dctI(evpt)*dctfactor;
figure
plot(t, evmiut)
hold on
plot(t, evpt, 'r')
max(abs(imag(evpt)))
figure
plot(w(1:1001), evmiuw(1:1001))
hold on
plot(w(1:1001), evpw(1:1001), 'r')
plot(w(1:1001), evpw(1:1001)./w(1:1001), 'r')
iw1000 = instwcos(psi(:, end));
iw1000 = instwcos(psi(:, end), 160);
figure
plot(x, iw1000)
iw1000f = instw(psi(:, end), 160);
figure
plot(x, iw1000f)
iw1000f = instw(psi(:, end), dx);
hold on
plot(x, iw1000f)
plot(x, iw1000f, 'r')
doc hilbert
iw1000f = instw(psi(:, end), dx);
figure
plot(0:254, phi)
plot(0:255, phi)
phi(2)-phi(1)
(phi(2)-phi(1))/pi
iw1000f = instw(psi(:, end), dx);
plot(0:255, phi)
phi = phi_cont(phi);
doc unwrap
doc round
mod(1.1, 1)
mod(-1.1, 1)
rem(-1.1, 1)
edit find
phi = angle(hilbert(fi0));
figure
plot(x, phi)
phi = phi_cont1(phi);
plot(x, phi)
phi = phi_cont1(phi);
plot(x, phi)
phi = angle(hilbert(fi0));
figure
plot(x, phi)
phic = phi_cont1(phi);
plot(x, phic)
phi = angle(hilbert(real(psi(:, end))));
figure
plot(x, phic)
plot(x, phi)
tic, phic = phi_cont1(phi); toc
tic, phic0 = phi_cont0(phi); toc
tic, phic0 = phi_cont(phi); toc
tic, phic = phi_cont1(phi); toc
figure
plot(x, phic)
plot(x, phic0)
iw1000 = instwcos(psi(:, end), 160);
figure
plot(x, iw1000)
iw1000f = instw(psi(:, end), dx);
hold on
plot(x, iw1000f, 'r')
edit hilbert
sign(1)
sign(-1)
sign(-9)
iw1000 = instwcos(psi(:, end), 160);
figure
plot(x, iw1000)
iw1000 = instwcos(psi(:, end), 160);
figure
plot(0:255, phi)
(phi(2)-phi(1))/pi
iw1000 = instwcos(psi(:, end), 160);
figure
plot(0:255, phi)
pi/2
(phi(2)-phi(1))/pi
(phi(2)-phi(1))
dphi(1)/pi*2
pi_half_jumps
(phi(206)-phi(205))/pi
iwi1000 = instwcos(imag(psi(:, end)), 160);
figure
plot(x, iwi1000)
iwi1000 = instwcos(imag(psi(:, end)), 160);
figure
plot(0:255, phi)
iwi1000 = instwcos(imag(psi(:, end)), 160);
plot(0:255, phi)
figure
plot(x, iwi1000)
iwi1000f = instw(imag(psi(:, end)), dx);
hold on
plot(x, iwi1000f, 'r')
%-- 19/05/2014 13:50 --%
rgb=imread('TSM2.jpg');
whos rgb
imshow(rgb)
truesize
gb=imread('TSM1.png');
rgb=imread('TSM1.png');
whos
imshow(rgb)
truesize
%-- 20/05/2014 13:05 --%
%-- 21/05/2014 14:42 --%
rgb=imread('HHGspectrum.jpeg');
whos rgb
imshow(rgb)
truesize
load coulomb_optV40
figure
plot(x, fi0)
plot(x(129:256), fi0(129:256))
figure
plot(x, fi0)
theta = 0:2*pi/128:(2*pi*127/128);
clear theta
data0 = real(fi0)*ones(1,128);
load tol4
data0 = real(fi0)*ones(1,256);
dataT = psi(:,end)*ones(1,256);
figure
mesh(data0)
mesh(dataT)
mesh(real(dataT))
whos
save twoDdata data0 dataT
clear all
load twoDdata
whos
%-- 26/05/2014 14:08 --%
load tol4
figure
plot(w(1:1001), evmiuw(1:1001))
whos
convtot = [conv, conv1(2:end)];
conv1(1)
conv(end)
figure
size(convtot)
convtot(end)-J
convtot(end)-J1
whos
figure
plot(0:184, convtot)
size(conv)
size(conv1)
%-- 29/05/2014 16:42 --%
x=0:2*pi/128:127*2*pi/128;
figure
y=-1:0.01:1;
plot(y, acos(y))
plot(y, asin(y))
r=asin((x/pi - 1));
figure
plot(x, r)
v=cos(x);
Dv = DfourierMap(v, r, 2*pi);
figure
plot(x, Dv)
Dv = DfourierMap(v, r, 2*pi);
figure;plot(r, cos(r))
r=asin((x/pi - 1))+pi/2;
plot(x, r)
r(end)+r(2)
plot(r, cos(r)
plot(r, cos(r))
Dv = DfourierMap(v, r, 2*pi);
plot(r, Dv)
v=cos(r);
Dv = DfourierMap(v, r, 2*pi);
plot(r, Dv)
r=(asin((x/pi - 1))+pi/2)*2;
v=cos(r);
plot(r, v)
Dv = DfourierMap(v, r, 2*pi);
plot(r, v)
plot(r, Dv)
theta=0:pi/128:127*pi/128;
r1=(-cos(theta) + 1)*2*pi;
figure
plot(r1, cos(r1))
r1=(-cos(theta) + 1)*pi;
plot(r1, cos(r1))
figure
plot(0, r1)
plot(r1, 0, '*')
v1=cos(r1);
Dv1 = DfourierMap(v1, r1, 2*pi);
figure
plot(r1, Dv1)
max(abs(Dv1+sin(r1)))
max(abs(Dv1+sin(r1.')))
max(Dv1)
Dv2 = DfourierMap(v1, r1, 2*pi);
figure
plot(r, Dv2)
Dv2 = DfourierMap(v1, r1, 2*pi);
plot(r, Dv2)
plot(r1, Dv2)
max(abs(Dv1+sin(r1.')))
max(abs(Dv1(2:127)+sin(r1(2:127).')))
max(abs(Dv2+sin(r1.')))
%-- 02/06/2014 11:05 --%
x=0:2*pi/128:127*2*pi/128;
r=rand(1, 128)*2*pi;
figure
plot(x, r)
r=sort(r);
plot(x, r)
plot(r, cos(r))
Dv = DfourierMap(cos(r), r, 2*pi);
figure
plot(r, Dv)
theta=0:pi/128:127*pi/128;
r1=(-cos(theta) + 1)*pi;
Dv = DfourierMap(cos(r1), r1, 2*pi);
figure
plot(r1, Dv)
max(abs(Dv+sin(r1.')))
max(abs(Dv(2:end-1)+sin(r1(2:end-1).')))
max(abs(Dv(3:end-2)+sin(r1(3:end-2).')))
hold on
plot(r1, sin(r1), 'r')
plot(r1, -sin(r1), 'r')
theta=0:pi/128:2*pi;
r1=(-cos(theta) + 1)*pi;
Dv = DfourierMap(cos(r1), r1, 2*pi);
figure
plot(r1, Dv)
theta=0:pi/128:pi;
r1=(-cos(theta) + 1)*pi;
Dv = DfourierMap(cos(r1), r1, 2*pi);
figure
plot(r1, Dv)
max(abs(Dv(2:end-1)+sin(r1(2:end-1).')))
max(abs(Dv+sin(r1.')))
Dv(1)
Dv = DfourierMap(cos(r1), r1);
r=(-cos(theta) + 1)*2*pi;
r=(-cos(x) + 1)*2*pi;
Dv = DfourierMap(cos(r), r);
plot(r, Dv)
hold on
plot(r, -sin(r), 'r')
max(abs(Dv+sin(r.')))
r=(-cos(x) + 1)*pi;
Dv = DfourierMap(cos(r), r);
clf
plot(r, Dv)
max(abs(Dv+sin(r.')))
plot(r, -sin(r), 'r')
max(abs(Dv(2:end-1)+sin(r(2:end-1).')))
max(Dv)
Dv(1)
Dv(2)
Dv(end)
Dv.'
figure
plot(r,Dv+sin(r), 'r')
plot(r,Dv+sin(r.'), 'r')
Dv = DcosinesMap(cos(r1), r1);
figure
plot(r1, Dv)
max(abs(Dv+sin(r1.')))
size(Dv)
size(r1)
max(abs(Dv+sin(r1)))
r2 = (theta/pi).^2;
figure
plot(theta, r2)
r2 = (theta/pi).^2*pi;
plot(theta, r2)
Dv = DfourierMap(cos(r2), r2);
Dv = DfourierMap(cos(r2(1:end-1)), r2(1:end-1));
plot(r2, Dv)
plot(r2(1:end-1), Dv)
r3=(1-cos(r2))*pi/2;
figure
plot(theta, r3)
Dv = DcosinesMap(cos(r3), r3);
figure
plot(r3, Dv)
hold on
plot(r3, -sin(r3), 'r')
max(abs(Dv+sin(r3)))
figure
plot(Dv+sin(r3))
plot(r3,Dv+sin(r3))
save
%-- 17/06/2014 10:32 --%
figure
whos
plot(t, fieldt)
doc plotyy
tfs = 2.41888e-17*t;
tfs(end)
tfs = 2.41888e-2*t;
tfs(end)
1000/40
ax1=gca
ax2=axes('Xaxislocation', 'top')
ax2=axes('position', get(ax1, 'position'),'Xaxislocation', 'top')
plot(tfs, fieldt)
clf
line(t, fieldt, 'color', 'k')
ax1=gca
ax2=axes('position', get(ax1, 'position'),'Xaxislocation', 'top', 'Yaxislocation', 'left')
ax2=axes('position', get(ax1, 'position'),'Xaxislocation', 'top', 'Yaxislocation', 'right')
figeps=gca
figeps('xdatasource', tfs)
tfs(end)
%-- 22/06/2014 12:58 --%
doc fsolve
%-- 23/06/2014 13:41 --%
syms a theta
int(a)
doc diff
int(a^2)
diff(a^2)
syms x y
diff(x*y, 2)
f(x)=x^2
syms x
f(x) = sin(x^2);
df = diff(f)
solve(x^2-1)
solve(x^2)
solve(x^2-1)
solve('x^2-1')
%-- 02/07/2014 12:16 --%
doc fsolve
tic, [r, Nr, Q] = map_grid1(@(r) tanh(r), 100); toc
tic, [r, Nr, Q] = map_grid1(@(r) tanh(r)+r, 100); toc
figure
plot(Q, r)
plot(r, Q)
Nr
tic, [r, Nr, Q] = map_grid1(@(r) tanh(r)+0.1*r, 100); toc
plot(r, Q)
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, 100); toc
plot(r, Q)
Nr
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, 100); toc
profiler
profile
profile(map_grid1(eqfun, maxr))
profile([r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, 100))
profile(map_grid1(@(r) 20*tanh(r) + r, 100))
[r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, 100)
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, 100); toc
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, @(r) 20*sech(r)^2 + 1 ,100); toc
figure
plot(r, Q)
tic, [r, Nr, Q] = map_grid1(@(r) tanh(r), @(r) sech(r)^2 ,100); toc
figure
plot(r, Q)
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r), @(r) 20*sech(r)^2 ,100); toc
plot(r, Q)
Nr
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r), @(r) 20*sech(r)^2 ,100); toc
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r) ,100); toc
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r), @(r) 20*sech(r)^2 ,100); toc
figure
plot(r, Q)
plot(Q, r)
plot(r, Q)
hold on
plot(r, 20*tanh(r), 'r')
tic, [r, Nr, Q] = map_grid1(@(r) 20*tanh(r) + r, @(r) 20*sech(r)^2 + 1 ,100); toc
clear all
%-- 06/07/2014 11:10 --%
x=-5:0.01:4.99;
v=exp(-x.^2/2);
figure
plot(x, v)
xnew=(rand(1, 128) - 1)*10;
hold on
vnewa = exp(-xnew.^2/2);
plot(xnew, vnewa)
xnew=sort(xnew);
vnewa = exp(-xnew.^2/2);
plot(xnew, vnewa, 'r')
xnew=(rand(1, 128) - 0.5)*10;
vnewa = exp(-xnew.^2/2);
plot(xnew, vnewa, 'r')
xnew=sort(xnew);
vnewa = exp(-xnew.^2/2);
plot(xnew, vnewa, 'r')
xdomain=[-5 5];
vnew = FourierIpln(v, xdomain, xnew);
max(abs(vnew-vnewa))
figure
plot(xnew, vnew-vnewa)
x2=0:2*pi/128:2*pi*127/128;
v2=cos(x2);
figure
plot(x2,v2)
xnew1=(rand(1, 128))*2*pi;
xnew1=sort(xnew1);
vnewa1 = cos(xnew1);
vnew1 = FourierIpln(v1, [0 2*pi], xnew1);
vnew1 = FourierIpln(v2, [0 2*pi], xnew1);
max(abs(vnew1-vnewa1))
x=-5:0.01:5;
v=exp(-x.^2/2);
vnew = cosineIpln(v, xdomain, xnew);
figure
max(abs(vnew1-vnewa1))
max(abs(vnew-vnewa))
plot(xnew, vnew)
plot(xnew, vnew-vnew(1))
max(abs(vnew-vnew(1)-vnewa))
vnew(1)
vnew = cosineIpln(v, xdomain, xnew);
max(abs(vnew-vnewa))
vnew1 = cosineIpln(v2, [0 2*pi], xnew1);
max(abs(vnew1-vnewa1))
figure
plot(xnew1, vnew1)
hold on
plot(xnew1, vnewa1, 'r')
vnew1 = cosineIpln(v2, [0 2*pi], xnew1);
x2=0:2*pi/128:2*pi;
v2=cos(x2);
vnew1 = cosineIpln(v2, [0 2*pi], xnew1);
max(abs(vnew1-vnewa1))
result = NewtonIpln4(x, v, xnew);
result = NewtonIpln4(x.', v.', xnew.');
tic, result = NewtonIpln4(x.', v.', xnew.'); toc
figure
plot(xnew, result)
size(result)
result = NewtonIpln4(x.', v.', xnew.');
result = NewtonIpln4b(x.', v.', xnew.');
size(result)
result = NewtonIpln4b(x.', v.', xnew.');
result = NewtonIpln4b(x, v, xnew.');
result = NewtonIpln4b(x, v, xnew);
size(result)
plot(xnew, result)
result = NewtonIpln4b(x(1:20:end), v(1:20:end), xnew);
plot(xnew, result)
result = NewtonIpln4b(x(1:50:end), v(1:50:end), xnew);
plot(xnew, result)
size(x)
result = NewtonIpln4b(x(1:100:end), v(1:100:end), xnew);
plot(xnew, result)
result = NewtonIpln4b(x(1:200:end), v(1:200:end), xnew);
plot(xnew, result)
%-- 08/07/2014 18:21 --%
x=-5:0.01:5;
v=exp(-x.^2/2);
result = NewtonIpln4(sp, fsp, 0.005);
result = NewtonIpln4(x, v, 0.005);
result
x2=0:2*pi/128:2*pi;
v2=cos(x2);
result = NewtonIpln4(x2, v2, pi+0.02);
result
cos(pi+0.02)
2*pi/128
result = NewtonIpln4(x2, v2, 0.02)
%-- 09/07/2014 16:01 --%
x=0:2*pi/6:2*pi;
v=cos(x);
result = NewtonIpln4(x, v, 2)
cos(2)
result = NewtonIpln4(x, v, 0.5)
cos(0.5)
result = NewtonIpln4b(x, v, (0:2*pi/128:2*pi).');
result = NewtonIpln4b(x.', v.', (0:2*pi/128:2*pi).');
result = NewtonIpln4b(x, v, (0:2*pi/128:2*pi));
size(result)
figure
x2=0:2*pi/128:2*pi;
plot(x2, result)
result = NewtonIpln4b(x, v, (0:2*pi/128:2*pi));
doc spdiags
result = NewtonIpln4b(x, v, (0:2*pi/128:2*pi));
size(result)
plot(x2, result)
hold on
plot(x2, cos(x2), 'r')
figure
plot(x2, result - cos(x2))
result = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
v2=cos(x2);
result = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
figure
plot(0:0.01:2*pi, result)
doc sort
[u ui]=sort([3 2 7 1 9])
u(ui)
uio=1:5;
u(uio(ui))
uir = uio(ui)
uir(ui) = uio
u(uir)
rand_order(1:10)
[u, ui] =rand_order(1:10)
u(ui)
[u, ui] =rand_order(1:10)
u(ui)
%-- 10/07/2014 12:03 --%
x2=0:2*pi/128:2*pi;
v2=cos(x2);
result = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
figure
plot(x2, result - cos(x2))
size(result)
plot(0:0.01:2*pi, result)
result = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
figure
plot(0:0.01:2*pi, result)
result = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
plot(0:0.01:2*pi, result)
figure
x3= 0:0.01:2*pi;
plot(x3, result-cos(x3))
resulto = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
figure
plot(x3, result-cos(x3))
plot(x3, resulto-cos(x3))
result = NewtonIpln4b(x2, v2, 0:0.01:2*pi);
figure
plot(x3, result-cos(x3))
x=-5:0.01:5;
v=exp(-x.^2/2);
vnew = NewtonIpln4b(x, v, -5:0.003:5);
figure
xnew = -5:0.003:5;
plot(xnew, vnew)
vnew = NewtonIpln4b(x, v, -5:0.003:5);
plot(xnew, vnew)
vnew = NewtonIpln4b(x(1:10:end), v, -5:0.003:5);
vnew = NewtonIpln4b(x(1:10:end), v(1:10:end), -5:0.003:5);
plot(xnew, vnew)
plot(xnew, vnew-exp(-x(1:10:end).^2/2))
plot(xnew, vnew-exp(-xnew.^2/2))
clear all
doc trapz
clear all
load('tol4.mat')
whos
fun = rectanglefun(w, 0.53, 0.55);
figure
plot(w, fun)
fun = rectanglefun(w, 0.53, 0.55).*evmiuw./(rectanglefun(w, 0.53, 0.55).*evmiuw.^2 + 1e-10);
figure
plot(w, fun)
plot(w(1:1001), fun(1:1001))
dw
dw=pi/1000
0.53/dw
0.55/dw
plot(w(160:180), fun(160:180))
fun = exp(-(w-0.54).^2/2*0.01^2);
plot(w(160:180), fun(160:180))
fun = exp(-(w-0.54).^2/(2*0.01^2));
plot(w(160:180), fun(160:180))
fun = exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw./(exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw.^2 + 1e-10);
plot(w(160:180), fun(160:180))
plot(w(100:300), fun(100:300))
figure
plot(t, dctI(fun))
fun = rectanglefun(w, 0.53, 0.55).*evmiuw./(rectanglefun(w, 0.53, 0.55).*evmiuw.^2 + 1e-10);
plot(t, dctI(fun))
fun = exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw./(exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw.^2 + 1e-10);
plot(t, dctI(fun))
fun = exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw./(exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw.^2 + 1e-8);
plot(w(100:300), fun(100:300))
fun = exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw./(exp(-(w-0.54).^2/(2*0.01^2)).*evmiuw.^2 + 1e-5);
plot(w(100:300), fun(100:300))
plot(t, dctI(fun))
dctfactor
figure
plot(t, dctI(rectanglefun(w, 0.53, 0.55).*evmiuw))
clear all
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
load coulomb_optV40
whos
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), evmiu(1:1001))
plot(w(1:1001), evmiuw(1:1001))
figure
plot(w(1:101), fieldw(1:101))
t=0:0.2:1e3;
plot(t, fieldt)
plot(t, evmiut)
plot(w(1:1001), log(evmiuw.^2(1:1001)))
plot(w(1:1001), log(evmiuw(1:1001).^2))
%-- 13/07/2014 12:52 --%
1300*12/365.25
ans*12
1300*12/365.25*13
t=0:0.1:100;
dw=pi/100;
w=0:dw:pi/0.1;
fw=exp(-(w-3).^2);
figure
plot(w, fw)
fieldwunc=zeros(1, 1001);
fieldwunc(96) = 1;
dctfactor = 100/sqrt(1e3*pi)
fieldtunc=dctI(fieldwunc)/dctfactor;
figure
plot(t, fieldtunc)
Nt = 1000;
vfilterE=fw;
dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilterE(1), vfilterE(2:Nt), 0.5*vfilterE(Nt + 1)])*dw;
%     dctfilterE = dctIintgrid([0.5*vfilterE(1); vfilterE(2:Nt); 0.5*vfilterE(Nt + 1)], T, t_ts(1:(Nt_ts-1)))/dctfactor;
coswT = ones(1, Nt + 1);
coswT(2:2:(Nt + 1)) = -1;
dctfilterET = sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1), vfilterE(2:Nt).*coswT(2:Nt), 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw;
deTdctfilterE  = dctfilterE0^2 - dctfilterET(allt_lasti)^2;
deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
lambda0 = (dctfilterE0*fieldtunc(1) - dctfilterET*fieldtunc(end))/deTdctfilterE;
lambdaT = (dctfilterE0*fieldtunc(end) - dctfilterET*fieldtunc(1))/deTdctfilterE
lambda0
fieldw = fieldwunc - fw.*(lambda0 + lambdaT*coswT);
figure
plot(w, fieldw)
fieldt = dctI(fieldw)/dctfactor;
plot(t, fieldt)
hold on
plot(t, fieldtunc, 'r')
figure
plot(t, dctI(fw)/dctfactor)
plot(t, dctI(fw).*coswT/dctfactor)
plot(t, dctI(fw.*coswT)/dctfactor)
fw2=rectanglefun(w, 2, 4);
figure
plot(w, fw2)
plot(t, dctI(fw2)/dctfactor)
plot(t, dctI(fw2.*coswT)/dctfactor)
dctfilterE02 = sqrt(2/pi)*sum([0.5*fw2(1), fw2(2:Nt), 0.5*fw2(Nt + 1)])*dw;
dctfilterET2 = sqrt(2/pi)*sum([0.5*fw2(1)*coswT(1), fw2(2:Nt).*coswT(2:Nt), 0.5*fw2(Nt + 1)*coswT(Nt + 1)])*dw;
deTdctfilterE2  = dctfilterE02^2 - dctfilterET2^2;
lambda02 = (dctfilterE02*fieldtunc(1) - dctfilterET2*fieldtunc(end))/deTdctfilterE2;
lambdaT2 = (dctfilterE02*fieldtunc(end) - dctfilterET2*fieldtunc(1))/deTdctfilterE2
fieldw2 = fieldwunc - fw2.*(lambda02 + lambdaT2*coswT);
figure
plot(w, fieldw2)
fieldt2 = dctI(fieldw2)/dctfactor;
plot(t, fieldt2)
hold on
plot(t, fieldtunc, 'r')
figure
plot(t, fieldtunc-fieldt)
plot(t, fieldtunc-fieldt2)
figure
plot(t, fieldt)
load toda
Vf = @(x) exp(-x) + x - 1;
[xdomain, dx, x] = nlVgrid(Vf, 128);
xdomain
[fi0, E0, x, E, P, H] = gsV(Vf, xdomain, 128);
E(1:10)
save toda xdomain dx fi0 E0 x E P H
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OClimfE0b(fi0, P(:, 5), Vf, [-10 85], xdomain, @(w) sech(20*(w-1).^4), @(w) 100*sech(20*(w-1).^4), 0.1, 100, 0.05, 5, 9, 1e-3);
figure
plot(0:niter, conv)
plot(t, fieldt)
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OClimfE0b(fi0, P(:, 5), Vf, [-10 85], xdomain, @(w) sech(20*(w-1).^4), @(w) 100*sech(20*(w-1).^4), 0.1, 100, 0.05, 5, 9, 1e-3);
figure
plot(0:0.05:100, dctI(fieldw)/dctfactor)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OClimfE0b(fi0, P(:, 5), Vf, [-10 85], xdomain, @(w) sech(20*(w-1).^4), @(w) 100*sech(20*(w-1).^4), 0.1, 100, 0.05, 5, 9, 1e-3);
w=0:dw:pi/0.05;
figure
plot(w, sech(20*(w-1).^4))
plot(w(1:101), sech(20*(w(1:101)-1).^4))
plot(w(1:101), sech(20*(w(1:101)-1).^4).*tanh(w(1:101)-1))
plot(w(1:101), sech(20*(w(1:101)-1).^4).*tanh(5*w(1:101)-1))
plot(w(1:101), sech(20*(w(1:101)-1).^4).*tanh(5*(w(1:101)-1)))
plot(w(1:101), sech(20*(w(1:101)-1).^4).*tanh(10*(w(1:101)-1)))
fw3=sech(20*(w(1:101)-1).^4);
figure
plot(w, dctI(fw3))
dctfactor = 100/sqrt(2e3*pi)
fw3=sech(20*(w-1).^4);
plot(w, dctI(fw3)/dctfactor)
plot(w, dctI(fw3.*tanh(10*(w-1)))/dctfactor)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OClimfE0b(fi0, P(:, 5), Vf, [-10 85], xdomain, @(w) sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 100*sech(20*(w-1).^4), 0.1, 100, 0.05, 5, 9, 1e-3);
figure
plot(0:niter, conv)
figure
t=0:0.05:100;
plot(t, fieldt)
plot(w(1:101), fieldw(1:101))
J1
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimfE0b([1;0], [0;1], [1 4], [-1 6], miu, @(w) @(w) sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 100*sech(20*(w-1).^4), 100, 0.1, 5, 5, 1e-3);
miu=[0 1; 1 0]
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimfE0b([1;0], [0;1], [1 4], [-1 6], miu, @(w) @(w) sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 100*sech(20*(w-1).^4), 100, 0.1, 5, 5, 1e-3);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimfE0b([1;0], [0;1], [1 4], [-1 6], miu, @(w) sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 100*sech(20*(w-1).^4), 100, 0.1, 5, 5, 1e-3);
figure
plot(0:niter, conv)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimfE0b([1;0], [0;1], [1 4], [-1 6], miu, @(w) sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 100*sech(20*(w-1).^4), 100, 0.1, 5, 5, 1e-3);
figure
plot(w, dctI(fw3)/dctfactor)
plot(w, dctI(fw3.*tanh(10*(w-1)))/dctfactor)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimfE0b([1;0], [0;1], [1 4], [-1 6], miu, @(w) 10*sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 100*sech(20*(w-1).^4), 100, 0.1, 5, 5, 1e-3);
J1
figure
plot(0:niter, conv)
t=0:0.1:100;
mLLNITERC
mallniterc
plot(t, fieldt)
plot(w(1:101), fieldw(1:101))
psiE=P\psi;
figure
plot(t, conj(psi).*psi)
plot(t, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfMLSE0b([1;0], [1 4], [-1 6], miu,...
@(w) sech(20*(w-1).^4).*tanh(10*(w-1)), @(w) 20*sech(20*(w-1).^4), @(w) exp(-10*(w-3).^2), 0.5, 100, 0.1, 5, 5, 1e-3);
figure
plot(0:niter, conv)
figure
plot(t, fieldt)
plot(t, fieldw)
plot(w, fieldw)
w=0:dw:pi/0.1;
figure
plot(t, evmiut)
plot(t, conj(psi).*psi)
J1
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 20*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 0.6*0.5*(1-tanh(100*(w-0.07))), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
%-- 21/07/2014 14:23 --%
figure
x=-5:0.01:5;
plot(x, erf(x))
plot(x, erf(x + 1i*x))
clear all
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.1*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
load coulomb_optV40
whos
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
figure
mx = evx(U, x);
figure
t=0:0.2:1e3;
plot(t, mx)
nU = normU(U);
figure
plot(t, nU)
mxw = dctI(mx)/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
mxw = dctI(mx)/dctfactor;
figure
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), mxw(1:1001))
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mxc = evx(Uc, x);
nUc = normU(Uc);
mxwc = dctI(mxc)/dctfactor;
figure
plot(t, mxc)
hold on;
plot(t, 0.03*sin(0.06*t), 'r')
hold on
plot(t, 0.03*cos(0.06*t), 'r')
mean(mx)
mean(mxc)
figure
plot(t, nUc)
figure
plot(w(1:1001), mxwc(1:1001))
viewVPmiux(U, Vabs, xabs, 0.06*sin(t), x, 0.01)
viewVPmiux(U(1:20:1e3), Vabs, xabs, 0.06*sin(t(1:20:1e3)), x, 0.01)
viewVPmiux(U(1:20:1e3), Vabs, xabs, 0.06*sin(t(1:20:1e3)), x, 0.1)
viewVPmiux(U, Vabs, xabs, 0.06*sin(t), x, 0.1)
viewVPmiux(Uc, Vabs, xabs, 0.06*cos(t), x, 0.1)
viewVlogPmiux(Uc, Vabs, xabs, 0.06*cos(t), x, 0.1)
viewVlogPmiux(U, Vabs, xabs, 0.06*sin(t), x, 0.1)
pi/(2*0.06)
figure
ans/0.2
plot(x, conj(U(:,132)).*U(:,132))
mx(132)
U26=U(:,132);
U26'*H*U26
ans-E0
whos
mp26 = U26'*ifft(p.*fft(U26));
mp26
mp = evp(U, p);
figure
plot(t, mp)
mpc = evp(Uc, p);
hold on
plot(t, mxc/10, 'r')
figure
plot(t, mpc)
ma = evmiu(U, -x./(1+x.^2).^(3/2) + 0.03*sin(0.06*t))
ma = evmiu(U, -x./(1+x.^2).^(3/2)) + 0.03*sin(0.06*t).*nU;
figure
plot(t, ma)
figure
plot(t, 0.03*cos(0.06*t)*0.5*(1 + tanh(t-5)))
plot(t, 0.03*cos(0.06*t).*0.5*(1 + tanh(t-5)))
plot(t, 0.03*cos(0.06*t).*0.5.*(1 + tanh(t-5)))
plot(t, 0.03*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-30))))
plot(t, 0.03*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-50))))
Eenv = 0.03*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-50)));
[Ue mnitere matvecse] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-50)), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[Ue mnitere matvecse] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-50))), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mxe = evx(Ue, x);
figure
plot(t, mxe)
mxwe = dctI(mxe)/dctfactor;
figure
plot(w(1:1001), mxwe(1:1001))
mpe = evp(Ue, p);
figure
plot(t, mpe)
mae = evmiu(Ue, -x./(1+x.^2).^(3/2)) + 0.03*sin(0.06*t).*nUe;
nUe = normU(Ue);
figure
plot(t, nUe)
nUe(end)
nU(end)
mae = evmiu(Ue, -x./(1+x.^2).^(3/2)) + 0.03*sin(0.06*t).*nUe;
figure
plot(t, mae)
hold on
plot(t, Eenv, 'r')
plot(t, Eenv/0.03*0.1849, 'r')
eqx = Vequilibrium(@(x) -x./(1+x.^2).^(3/2), 0.03*sin(0.06*t));
figure
plot(t, eqx)
hold on
plot(t, mx, 'r')
figure
whos
plot(x, Vf(x))
figure
plot(x, xabs)
plot(x, xabs.*0.5.*(tanh(20*(x-1)) - tanh(20*(x+1))))
plot(x, 0.5.*(tanh(20*(x-1)) - tanh(20*(x+1))))
plot(x, 0.5.*(tanh(20*(x-1)) - tanh(20*(x+1))+ 2))
plot(x, xabs*0.5.*(tanh(20*(x-1)) - tanh(20*(x+1))+ 2))
plot(x, Vf(x))
plot(x, xabs*0.5.*(tanh(10*(x-2)) - tanh(10*(x+2))+ 2))
xabsc = xabs*0.5.*(tanh(10*(x-2)) - tanh(10*(x+2))+ 2);
[Uc1 mniterc1 matvecsc1] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabsc*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
figure
mxc1 = evx(Uc1, x);
plot(x, mxc1)
plot(t, mxc1)
hold on
plot(t, mxc, 'r')
hold on
plot(x, fi0, 'r')
plot(x, Vf(x)-xabs*0.03)
E0
figure
plot(x, xabs*0.5.*(tanh((x-15)) - tanh((x+15))+ 2))
xabs2 = xabs*0.5.*(tanh(10*(x-2)) - tanh(10*(x+2))+ 2);
[Uc2 mniterc2 matvecsc2] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs2*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mxc2 = evx(Uc2, x);
figure
plot(t, mxc2)
hold on
plot(t, mxc1, 'r')
xabs2 = xabs*0.5.*(tanh((x-15)) - tanh((x+15))+ 2));
xabs2 = xabs*0.5.*(tanh((x-15)) - tanh((x+15))+ 2);
[Uc2 mniterc2 matvecsc2] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs2*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
clf
mxc2 = evx(Uc2, x);
plot(t, mxc2)
whos
%-- 22/07/2014 13:47 --%
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
load coulomb_optV40
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mxc = evx(Uc, x);
mpc = evp(Uc, p);
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(t, mxc)
plot(t, mpc)
hold on
plot(t, 0.03*cos(0.06*t), 'r')
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mx = evx(U, x);
mp = evp(U, p);
figure
plot(t, mp)
hold on
plot(t, 0.03*cos(0.06*t), 'r')
figure
plot(t, mp)
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.05*p.').*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.05*p).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.05*x).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mx1 = evx(U1, x);
figure
plot(t, mx1)
mp1 = evp(U1, p);
figure
plot(t, mp1)
mxw1 = dctI(mx1)/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
mxw1 = dctI(mx1)/dctfactor;
figure
plot(w(1:1001), mxw(1:1001))
plot(w(1:1001), mxw1(1:1001))
viewVlogPmiux(U1, Vabs, xabs, 0.06*sin(t), x, 0.1)
viewVPmiux(U1, Vabs, xabs, 0.06*sin(t), x, 0.1)
[U2 mniter2 matvecs2] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.04*x).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mx2 = evx(U2, x);
mxw2 = dctI(mx2)/dctfactor;
figure
plot(t, mx2)
[U2 mniter2 matvecs2] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.02*x).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
norm(U(:, end))
norm(U1(:, end))
norm(U2(:, end))
mx2 = evx(U2, x);
mxw2 = dctI(mx2)/dctfactor;
plot(t, mx2)
mx = evx(U, x);
mxw = dctI(mx)/dctfactor;
hold on
plot(w(1:1001), mxw(1:1001))
plot(w(1:1001), mxw(1:1001), 'r')
sinh(1)
k1=0:0.01:10;
D=4*k1.^2./((k1.^2-1).^2*sinh(1) + 4*k1.^2);
figure
plot(k1, D)
norm(Uc(:, end))
norm(U1(:, end))
D=4*k1.^2./((k1.^2+1).^2*sinh(1) + 4*k1.^2);
plot(k1, D)
E=k1.^2/2;
kapa=sqrt(2*(5-E));
kapa=sqrt(2*(10-E));
D=4*k1.^2.*kapa.^2./((k1.^2+kapa.^2).^2*sinh(kapa).^2 + 4*k1.^2.*kapa.^2);
D=4*k1.^2.*kapa.^2./((k1.^2+kapa.^2).^2.*sinh(kapa).^2 + 4*k1.^2.*kapa.^2);
plot(k1, D)
%-- 23/07/2014 12:32 --%
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.05*x).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
load coulomb_optV40
[U1 mniter1 matvecs1] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], exp(1i*0.05*x).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U10 mniter10 matvecs10] = TDHxpKr(K, Vabs, @(u,x,t) 0, [], exp(1i*0.05*x).*fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[Ue mnitere matvecse] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-50))), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mx = evx(U, x);
mxw = dctI(mx)/dctfactor;
mxe = evx(Ue, x);
mxc = evx(Uc, x);
mxwc = dctI(mxc)/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
mx1 = evx(U1, x);
mxw = dctI(mx)/dctfactor;
mxw1 = dctI(mx1)/dctfactor;
mxwc = dctI(mxc)/dctfactor;
mxwe = dctI(mxe)/dctfactor;
mxw10 = dctI(mx10)/dctfactor;
mx10 = evx(U10, x);
mxw10 = dctI(mx10)/dctfactor;
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(t, mx1)
hold on
plot(t, mx10, 'r')
figur
figure
plot(w(1:1001), mxw1(1:1001))
hold on
plot(w(1:1001), mxw10(1:1001), 'r')
norm(U10(:, end))
norm(U1(:, end))
UE1 = P\U1;
UE10 = P\U10;
figure
plot(t, UE1)
plot(t, conj(UE1).*UE1)
figure
plot(t, conj(UE10).*UE10)
UE = P\U;
figure
plot(t, conj(UE).*UE)
UEc = P\Uc;
plot(t, conj(UEc).*UEc)
UEe = P\Ue;
plot(t, conj(UEe).*UEe)
nUe = normU(Ue);
max(abs(nUe - sum(conj(Ue).*Ue)))
figure
plot(t, nUe - sum(conj(Ue).*Ue))
figure
plot(t, nUe)
nUe1=sum(conj(Ue).*Ue);
figure
plot(t, nUe1)
hold on
plot(t, nUe, 'r')
hold on
plot(t, nUe, 'r')
doc norm
max(abs(nUe - sqrt(sum(conj(Ue).*Ue))))
mE0 = evHt(U, K, Vabs);
mE = evHt(U, K, Vabs, 0.03*sin(0.06*t), xabs);
figure
plot(t, mE0)
figure
plot(t, mE)
max(mE)-min(mE)
max(mE0)-min(mE0)
hold on
plot(t, 0.03*sin(0.06*t), 'r')
plot(t, 0.03*sin(0.06*t)/10 + 0.3, 'r')
plot(t, 0.03*sin(0.06*t)/15 + 0.33, 'r')
mEc0 = evHt(Uc, K, Vabs);
figure
plot(t, mEc0)
hold on
plot(t, mE, 'r')
plot(t, mE0, 'r')
sqnc = sqnorm(Uc);
sqnc = sqnormpsi(Uc);
figure
plot(t, sqnc)
hold on
plot(t, 0.03*cos(0.06*t)/15 + 0.33, 'g')
mEe0 = evHt(Ue, K, Vabs);
figure
plot(t, mEe0)
mE10 = evHt(U1, K, Vabs);
figure
plot(t, mE10)
hold on
plot(t, mE0, 'r')
figure
plot(t, mxc)
2*pi/(37.4-20.6)
2*pi/(20.6-4.6)
figure
plot(t, mx1)
plot(t, sqnc)
2*pi/0.06
ans/4
Eenv = 0.03*cos(0.06*t).*0.5.*(1 + tanh(0.1*(t-50)));
733.2/0.2
[fi0733, E0733] = gsV(Vabs - x*Eenv(3667), xdomain, Nx);
E0733
figure
plot(t, mxe)
[fi0733, E0733] = gsV(Vabs - xabs*Eenv(3667), xdomain, Nx);
E0733
figure
plot(x, Vabs - xabs*Eenv(3667))
hold on
plot(x, fi0733)
whos
[fi0733, E0733] = gsV(Vf(x) - x*Eenv(3667), xdomain, Nx);
E0733
plot(x, fi0733)
[fi0733, E0733, ~, E733, P733] = gsV(Vf(x) - x*Eenv(3667), xdomain, Nx);
E733
[(1:15).' E733]
[(1:15).' E733(1:15)]
[(1:20).' E733(1:20)]
figure
plot(t, mxe)
plot(t, mEe)
mEe = evHt(Ue, K, Vabs, Eenv, xabs);
plot(t, mEe)
[(1:30).' E733(1:30)]
figure
plot(x, P(:, 24))
plot(x, P733(:, 24))
plot(x, P733(:, 24).*conj(P(:, 24)))
plot(x, P733(:, 24).*conj(P733(:, 24)))
plot(x, P733(:, 23).*conj(P733(:, 23)))
plot(x, P733(:, 25).*conj(P733(:, 25)))
plot(x, P733(:, 22).*conj(P733(:, 22)))
plot(x, P733(:, 26).*conj(P733(:, 26)))
plot(x, P733(:, 24).*conj(P733(:, 24)))
hold on
plot(x, Ue(:, 3667).*conj(Ue(:, 3667)), 'r')
max(abs(Ue(:, 3667).*conj(Ue(:, 3667))- P733(:, 24).*conj(P733(:, 24))))
figure
plot(x, Ue(:, 3667).*conj(Ue(:, 3667))- P733(:, 24).*conj(P733(:, 24)))
figure
plot(t, mEe)
916.2-864
2*pi/0.06
ans/2
2*pi/0.06/4
916.4+ans
2*pi/0.06/8
916.4+ans
hold on
plot(t, Eenv/15 + 0.33, 'g')
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalx(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-6, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxE0b(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-6, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, [-2, 15.5], xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-6, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-6, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 0.6*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-5, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 1.2*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-5, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-5, 2e3, 0.2, 5, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-5, 2e3, 0.2, 9, 9, 5e-4);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(fi0, Vf, 1, xdomain, @(x) x, @(x) 0.5*(tanh(x+35)-tanh(x-35)), @(x) 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1), @(w) 2*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e5*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-5, 2e3, 0.2, 9, 9, 5e-4);
sum(J1fun)
figure
plot(w, evmiualw)
plot(0:pi/2e3:pi, evmiualw(1:2001))
figure
plot(w(1:101), fieldw(1:101)
plot(w(1:101), fieldw(1:101))
w=0:pi/2e3:pi/dt;
plot(0:pi/2e3:pi/10, fieldw(1:201))
sum(J1fun)
figure
plot(0:0.2:2e3, fieldt)
plot(0:pi/2e3:pi/10, fieldw(1:201))
plot(0:pi/2e3:pi, evmiuw(1:2001))
sqnorm(psi(:, end))
J1
load HCl1
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1,5e-2), @(w) 10000*rectanglefun(w, 0, 1.5e-2),...
@(w) rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 10000*rectanglefun(w, 0, 1.5e-2),...
@(w) rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 10000*rectanglefun(w, 0, 1.5e-2),...
@(w) rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 10000*rectanglefun(w, 0, 1.5e-2),...
@(w) rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
figure
plot(0:5:1e4, fieldt)
figure
plot(0:pi/1e4:pi/100, fieldw(1:101))
E(2)-E(1)
figure
plot(0:pi/1e4:pi/100, evmiuw(1:101))
viewVPmiux(psi, Vf, x, fieldt, x, 0.1)
J1
hold on
plot(0:pi/1e4:pi/100, fieldw(1:101), 'r')
figure
plot(0:5:1e4, fieldt)
plot(0:5:1e4, fieldt, 'r')
hold on
plot(0:5:1e4, evmiut)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2),...
@(w) 100*rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2),...
@(w) 100*rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2),...
@(w) 100*rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1e-4, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiuxE0b(fi0, Vf, 1785, [-0.08 0.25],...
xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2)*cos(2*pi*w/1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2),...
@(w) 100*rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1e-3, 10000, 5, 5, 9, 1e-3);
J1
figure
plot(0:5:1e4, evmiut)
hold on
plot(0:5:1e4, fieldt, 'r')
figure
plot(w(1:101), fieldw(1:101), 'r')
hold on
plot(w(1:101), evmiuw(1:101))
viewVPmiux(psi, Vf, x, fieldt, x, 0.1)
save HClE0b fieldt fieldw evmiut evmiuw relE conv niter J1
iwfield = instwcos(fieldt, 1e4);
figure
plot(t, iwfield)
plot(0:5:1e4, iwfield)
iwfield1 = instw(fieldt, 5);
figure
plot(0:5:1e4, iwfield1)
hold on
plot(0:5:1e4, iwfield, 'r')
iwevmiu = instw(evmiut, 5);
plot(0:5:1e4, iwevmiu)
E(2:3)-E0
hold on
plot(0:5:1e4, evmiut, 'r')
psiE= P\psi;
figure
plot(0:5:1e4, conj(psiE).*psiE)
1515-1045
2*pi/ans
2465-1985
2*pi/ans
920-440
1330-440
1330-920
2*pi/ans
1780-1330
%-- 27/07/2014 14:49 --%
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
load coulomb_optV40
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sum(J1fun)
penalnormf(normpsiT)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sum(J1fun)
penalnormf(normpsiT)
sum(J1fun)
penalnormf(normpsiT)
sum(J1fun)
sum(J2fun)
J1
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
viewVPmiux(psi, Vf, xabs, fieldt, x, 0.1)
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
J, J1, J2, Jorth, Jpnorm
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 2*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
J, J1, J2, Jorth, Jpnorm
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
figure
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), evmiuw(1:1001))
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
figure
t=0:0.2:1e3;
plot(t, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 25*sech(50*(x-0.9))^2, @(w) 2*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
figure
psi(:,end)'*psi(:,end)
plot(t, fieldt)
figure
plot(w(1:1001), evmiuw(1:1001))
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07))/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07))/dctfactor;
hold on
plot(t, fieldtg, 'r')
max(abs(fieldt-fieldtg)
max(abs(fieldt-fieldtg))
figure
plot(fieldt - fieldtg)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 0.5e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 0.5e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 2e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 2e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldt)
hold on
plot(t, fieldtg, 'r')
J1
figure
plot(w(1:1001), evmiuw(1:1001))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
figure
plot(t, fieldtg, 'r')
hold on
plot(t, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 0.5*(1-tanh(100*(w-0.07))).*cos(2*pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J
]
[J, J1, J2, Jorth, Jpnorm]
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(w*1e3), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
figure
plot(t, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
figure
plot(t, fieldt)
[J, J1, J2, Jorth, Jpnorm]
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
figure
plot(t, fieldt)
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07))/dctfactor;
hold on
plot(t, fieldtg, 'r')
figure
plot(w(1:1001), evmiuw(1:1001))
whos
xdomain
clear all
t=0:0.2:1e3;
whos
w=0:pi/1e3:pi/0.2;
dctfactor = 1e3/(sqrt(5e3*pi))
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
load coulomb_optV120
whos
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0160, Vabs160, 1, xdomain, xabs160, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
xdomain=[-160 160]
x(1)
x160(1)
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0160, Vabs160, 1, xdomain, xabs160, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0160, Vabs160, 1, xdomain, xabs160, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0160, Vabs160, 1, xdomain, xabs160, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
psi(:,end)'*psi(:,end)
figure
plot(w(1:1001), evmiuw(1:1001))
figure
plot(t, fieldt)
%-- 29/07/2014 12:04 --%
load coulomb_optV40
whos
E(1:10)-E0
E(2:10)-E(1:9)
[fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.17, 0.19), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), evmiuw(1:1001))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.17, 0.19), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:dw:pi, fieldw(1:101))
plot(0:dw:pi, fieldw(1:1001))
plot(0:dw:pi/10, fieldw(1:101))
plot(0:dw:pi, evmiuw(1:1001))
plot(0:dw:pi/10, fieldw(1:101))
plot(0:0.2:1e3, dctI(fieldw)/dctfactor)
figure
plot(t, fieldt)
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07))/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07))/dctfactor;
hold on
plot(t, fieldtg, 'r')
figure
plot(w(1:101), fieldw(1:101))
hold on
fieldwg = 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07);
plot(w(1:101), fieldwg(1:101), 'r')
figure
plot(w(1:1001), evmiuw(1:1001))
psiE= P\psi;
figure
plot(t, psiE)
plot(t, conj(psiE).*psiE)
psi(:,end)'*psi(:,end)
J1
J
conv(end)
figure
plot(t, evmiut)
evxt = evmiu(psi, x);
hold on
plot(t, evxt, 'r')
evpt = evp(psi, p);
evwp = dctI(evpt)*dctfactor;
hold on
plot(w(1:1001), evpw(1:1001))
evpw=evwp;
clear evwp
plot(w(1:1001), evpw(1:1001), 'r')
figure
plot(w(1:1001), evpw(1:1001), 'r')
evat = eva(psi, -x./(1+x.^2).^(3/2), fieldt);
figure
plot(t, evpt)
plot(t, evat)
evaw = dctI(evat)*dctfactor;
figure
plot(w(1:1001), evaw(1:1001), 'r')
%-- 04/08/2014 14:30 --%
load coulomb_optV40
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -0.5*xabs*cos(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniterc
[Uc mniterc matvecsc] = TDHxpKr(K, Vabs, @(u,x,t) -xabs*cos(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
solveOCkr
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.03*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.3*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 5e3, 9, 9, 1e-5);
mniter
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.3*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
mniter
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.3*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 7, 7, 1e-5);
mniter
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.3*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 7, 7, 1e-5);
[U mniter matvecs] = TDHxpKr(K, Vabs, @(u,x,t) -0.3*xabs*sin(0.06*t), [], fi0, x, [0 1e3], 1e4, 9, 9, 1e-5);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnormE0b(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(t, fieldt)
hold on
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07))/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07))/dctfactor;
hold on
plot(t, fieldtg , 'r')
figure
plot(w(1:1001), evmiuw(1:1001))
fieldwg = 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07);
figure
plot(w(1:101), fieldwg(1:101), 'r')
hold on
plot(w(1:101), fieldw(1:101))
psiE= P\psiE;
psiE= P\psi;
[fieldtg, fieldwg, psig, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldtg, fieldwg, psig, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresultspnorm(fi0, Vabs, 1, xdomain, xabs, @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
hold on
plot(w(1:1001), evmiuwg(1:1001), 'r')
figure
plot(t, conj(psiE).*psiE)
psi(:,end)'*psi(:,end)
test_harc
testTLSmat
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
size(x)
size(fi0)
size(xabs)
clear all
load coulomb_optV40
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
M=spdiags(allfield, 0, allt_lasti, allt_lasti);
spdiags(allfield, 0, allt_lasti, allt_lasti);
diag(allfield);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
dctfactor = 1e3/(sqrt(5e3*pi))
figure
plot(t, evat)
figure
plot(t, fieldt)
figure
plot(w(1:1001), evmiuw(1:1001))
figure
plot(w(1:1001), evaw(1:1001))
hold on
plot(w(1:1001), evaw(1:1001)./w(1:1001).^2, 'r')
%-- 06/08/2014 19:59 --%
[y,k]=f(188,9)
[y,k]=f(188*0.05,9)
test
%-- 07/08/2014 07:56 --%
[y,k]=f(188*0.05,9)
[y,k]=f(188*0.05,12)
[y,k]=f(188*0.05,15)
[y,k]=f(188*0.05,20)
[y,k]=f(188*0.05,9)
[y,k]=f(-1i*188*0.05,9)
[y,k]=f(1i*188*0.05,9)
[y,k]=f(-1i*188*0.05,9)
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
load
clear all
load coulomb_optV40
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.1*25*sech(50*(x-0.9))^2, @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J1
[y,k]=f(-1i*188*0.05,9)
abs(y)
abs(-1i*188*0.05,9)
abs(-1i*188*0.05)
result = my_direct_f(-1i*188*0.05, 9)
y
[y,k]=f(100,9)
result = my_direct_f(100, 9)
result -y
result - y
y
y-result
result
y
(y-result)./result
[y,k]=f(10,9)
result = my_direct_f(10, 9)
result = my_direct_f(1i*10, 9)
[y,k]=f(1i*10,9)
[y,k]=f(1i*30,9)
result = my_direct_f(1i*30, 9)
[y,k]=f(1i*50,9)
result = my_direct_f(1i*50, 9)
result = my_direct_f(1i*100, 9)
[y,k]=f(1i*100,9)
result = my_direct_f(30, 9)
[y,k]=f(30,9)
[y,k]=f(0.05*188,9)
result = my_direct_f(0.05*188, 9)
result = my_direct_f(20, 9)
[y,k]=f(20,9)
[y,k]=f(1i*20,9)
result = my_direct_f(1i*20, 9)
result = my_direct_f(-1i*20, 9)
[y,k]=f(-1i*20,9)
[y,k]=f(-1i*188*0.05,9)
result = my_direct_f(-1i*188*0.05, 9)
(y-result)./result
f_comparison(9);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = geussresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jorth, Jpnorm]
[Jg, J1g, J2g, Jorthg, Jpnormg]
figure
w=0:pi/1e3:pi/0.2;
t=0:0.2:1e3;
plot(w(1:1001), evaw(1:1001))
plot(t, fieldtg)
plot(t, evag)
plot(t, evatg)
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 4*0.5*(1-tanh(100*(w-0.07))).*cos(pi*w/0.07), @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 1, 2), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
mniterc
mnitercg
%-- 07/08/2014 15:21 --%
bench
%-- 11/08/2014 07:49 --%
load coulomb_optV40
whos
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, weight1] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.38, 0.40), 0, weight, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), evaw(1:1001))
psi(:,end)'*psi(:,end)
psi1(:,end)'*psi1(:,end)
plot(w(1:1001), evaw1(1:1001))
psiE= P\psi;
psiE1= P\psi1;
figure
plot(t, conj(psiE1).*psiE1)
figure
plot(t, fieldt1)
figure
plot(t, evmiut1)
plot(t, eva1)
plot(t, evat1)
J11
J1
conv(end)
viewVPmiux(psi, Vf, xabs, fieldt1, x, 0.1)
viewVPmiux(psi1, Vf, xabs, fieldt1, x, 0.1)
viewVPmiux(psi1, Vf, xabs, fieldt1, x, 0.01)
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
mniterc
[J, J1, J2, Jpnorm]
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1, @(w) 1e4*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[J, J1, J2, Jpnorm]
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.61, 0.63), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.61, 0.63), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(w(1:1001), evaw(1:1001))
figure
plot(w(1:1001), evmiuw(1:1001))
psiE= P\psi;
figure
plot(t, conj(psiE).*psiE)
E(5)-E(1)
E(6)-E(1)
E(4)-E(1)
figure
E(13)-E(1)
E(13)-E(1:10)
plot(t, conj(psiE(6,:)).*psiE(6,:))
plot(t, conj(psiE(5,:)).*psiE(5,:))
plot(t, conj(psiE(3:7,:)).*psiE(3:7,:))
plot(t, conj(psiE(2:7,:)).*psiE(2:7,:))
figure
plot(t, fieldt1)
hold on
plot(t, fieldt, 'r')
J1
hold on
plot(w(1:1001), evmiuw1(1:1001), 'r')
plot(w(1:1001), evaw1(1:1001), 'r')
figure
plot(t, fieldt1 - fieldt)
plot(t, (fieldt1 - fieldt./fieldt))
plot(t, (fieldt1 - fieldt)./fieldt)
plot(t, (fieldt1 - fieldt)./(fieldt+eps))
plot(t, (fieldt1 - fieldt)./(fieldt+1e-5))
plot(t, (fieldt1 - fieldt)./(fieldt+1e-4))
plot(t, (fieldt1 - fieldt)./(fieldt+1e-3))
plot(t, (fieldt1 - fieldt)./(fieldt+1e-2))
figure
plot(t, conj(psiE(2:7,:)).*psiE(2:7,:))
figure
plot(t, evat)
plot(t, evmiut)
hold on
plot(t, evmiut1, 'r')
figure
plot(w(1:1001), evaw1(1:1001))
plot(w(1:1001), emiuw1(1:1001))
plot(w(1:1001), evmiuw1(1:1001))
figure
plot(w(1:1001), evaw1(1:1001))
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, weight2] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.8, 1.2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) rectanglefun(w, 0.8, 1.2), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, weight2] = OCfpnorm_evaE0b(fi0, Vabs, 1, xdomain, xabs, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.8, 1.2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, field1)
plot(t, fieldt1)
hold on
plot(t, fieldt2, 'r')
figure
plot(w(1:1001), evaw2(1:1001))
psi(:,end)'*psi(:,end)
figure
plot(t, evat2)
plot(t, evmiut2)
viewVPmiux(psi2, Vf, xabs, fieldt2, x, 0.01)
psiE2= P\psi2;
figure
plot(t, conj(psiE2).*psiE2)
figure
plot(t, evat2)
figure
psi1(:,end)'*psi1(:,end)
psi2(:,end)'*psi2(:,end)
plot(t, conj(psiE1).*psiE1)
figure
plot(t, log10(evat2))
plot(w(1:2001), log10(evat2(1:2001))
plot(w(1:2001), log10(evat2(1:2001)))
plot(w, log10(evat2))
plot(w, log10(abs(evat2)))
plot(w, log10(evat2.^2))
plot(w, log10(evat2))
plot(w, log10(abs(evat2)))
plot(w, log10(abs(evaw2)))
whos
[fieldt160, fieldw160, psi160, evat160, evaw160, evmiut160, evmiuw160, relE160, conv160, niter160, mallniterc160, J1160, maxgrad160, weight160] = OCfpnorm_evaE0b(fi0160, Vabs160, 1, [-160 160], xabs160, -x./(1 + x.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.8, 1.2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
xdomain
xdomain160=[-160 160];
[fieldt160, fieldw160, psi160, evat160, evaw160, evmiut160, evmiuw160, relE160, conv160, niter160, mallniterc160, J1160, maxgrad160, weight160] = OCfpnorm_evaE0b(fi0160, Vabs160, 1, [-160 160], xabs160, -x160./(1 + x160.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.8, 1.2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
psi160(:end)'*psi160(:,end)
psi160(:,end)'*psi160(:,end)
plot(t, fieldt1)
hold on
plot(t, fieldt2, 'r')
plot(t, fieldt160, 'g')
figure
plot(w(1:1001), evaw160(1:1001))
viewVPmiux(psi160, Vf, xabs160, fieldt160, x, 0.01)
viewVPmiux(psi160, Vf, xabs160, fieldt160, x160, 0.01)
figure
plot(t, evat160)
figure
plot(t, evmiut160)
sqn160 = sqnorm(psi160);
figure
plot(t, sqn160)
whos
psiE160=P160\psi160;
figure
plot(t, conj(psiE160).*psiE160)
save fieldw_fundamental fieldw1
load coulomb_optV240
[fieldt240, fieldw240, psi240, evat240, evaw240, evmiut240, evmiuw240, relE240, conv240, niter240, mallniterc240, J1240, maxgrad240, weight240] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 20*rectanglefun(w, 0.8, 1.2), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
save teu fieldw1 fieldw fieldw160 fieldw240
figure
plot(w(1:1001), evaw160(1:1001))
hold on
plot(w(1:1001), evaw240(1:1001))
plot(w(1:1001), evaw240(1:1001),'r')
psi160(:,end)'*psi160(:,end)
psi240(:,end)'*psi240(:,end)
J1160
J1240
figure
plot(t, evmiut160)
hold on
plot(t, evmiut240, 'r')
figure
plot(t, evat160)
hold on
plot(t, evat240, 'r')
figure
plot(t, fieldt160)
hold on
plot(t, fieldt240, 'r')
%-- 14/08/2014 17:37 --%
whos
load coulomb_optV240
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw31, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 1e4*rectanglefun(w, 1.85, 1.87), 0, 1e3, 0.2, 7, 7, 1e-3);
Jg
conv31(end)
mnitercg
[Jg, J1g, J2g, Jorthg, Jpnormg]
sqnorm(psig(:,end))
figure
plot(w(1:1001), evaw31(1:1001))
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), evaw31(1:1001))
plot(w(1:1001), evawg(1:1001))
figure
plot(t, fielt31)
plot(t, fieldt31)
plot(t, fieldtg)
hold on
dctfactor = 1e3/(sqrt(5e3*pi))
fieldt2=dctI(fieldw)/dctfactor;
fieldt2=dctI(fieldw2)/dctfactor;
plot(t, fieldt2, 'r')
[fieldt31, fieldw31, psi31, evat31, evaw31, evmiut31, evmiuw31, relE31, conv31_2, niter31_2, mallniterc31, J131, maxgrad31, weight31] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw31, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 1e4*rectanglefun(w, 1.85, 1.87), 0, weight31, 1e3, 0.2, 7, 7, 1e-3, 1e3);
J131
figure
plot(w(1:1001), evaw31(1:1001))
[fieldt2, fieldw2_2, psi2, evat2, evaw2, evmiut2, evmiuw2, mniterc2, J2, J12, J22, Jorth2, Jpnorm2] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw2, @(w) 1e3*0.5*(1-tanh(100*(w-0.07))), @(w) 2000*rectanglefun(w, 1.2, 2), 0, 1e3, 0.2, 7, 7, 1e-3);
J2
conv2(end)
J22
[J12
J22 Jpnorm2]
[J12 J22 Jpnorm2]
J2
J12/100+J22+Jpnorm
J12/100+J22+Jpnorm2
J12/10+J22+Jpnorm2
J12/10+J22*10+Jpnorm2
conv2(end)
[fieldt2, fieldw2_2, psi2, evat2, evaw2, evmiut2, evmiuw2, mniterc2, J2, J12, J22, Jorth2, Jpnorm2] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw2, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 200*rectanglefun(w, 1.2, 2), 0, 1e3, 0.2, 7, 7, 1e-3);
[J12 J22 Jpnorm2]
J2
figure
plot(w(1:1001), evaw2(1:1001))
41*0.06
31*0.06
35*0.06
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw2, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 1e3*rectanglefun(w, 2.09, 2.11), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw2, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 2.09, 2.11), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35_2, niter35_2, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw2, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 2.09, 2.11), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt35, fieldw35, psi35, evat35, evaw35, evmiut35, evmiuw35, relE35, conv35_2, niter35_2, mallniterc35, J135, maxgrad35, weight35] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw2, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 2.09, 2.11), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sqnorm(psig(:,end))
sqnorm(psi2(:,end))
figure
plot(t, fieldt2)
hold on
plot(t, fieldt35, 'r')
figure
plot(w(1:1001), evaw2(1:1001))
hold on
plot(w(1:1001), evaw35(1:1001), 'r')
figure
plot(w(1:2001), evaw35(1:2001))
figure
plot(t, evat2)
hold on
plot(t, eva35, 'r')
plot(t, evat35, 'r')
sqnorm(psi2(:,end))
sqnorm(psi35(:,end))
npsi2=sqnorm(psi2);
npsi35=sqnorm(psi35);
figure
plot(t, npsi2)
hold on
plot(t, npsi35, 'r')
viewVPmiux(psi35, Vf, xabs240, fieldt35, x240, 0.01)
viewVPmiux(psi35, Vabs240, xabs240, fieldt35, x240, 0.01)
figure
plot(t, evmiut31)
plot(t, evmiut35)
figure
plot(w(1:2001), evmiuw35(1:2001))
save field35 fieldt35 fieldw35 evat35 evaw35 evmiut35 evmiuw35 conv35_2 niter35_2 J135
%-- 15/08/2014 14:46 --%
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(w(1:101), exp(-(w-0.06)^2/(2*0.03^2)))
plot(w(1:101), exp(-(w(1:101)-0.06).^2/(2*0.03^2)))
plot(w(1:101), exp(-(w(1:101)-0.06).^2/(2*0.02^2)))
plot(w(1:101), exp(-(w(1:101)-0.06).^2/(2*0.01^2)))
load coulomb_optV240
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
figure
plot(w(1:101), fieldw)
plot(w(1:101), fieldwg)
plot(w(1:101), fieldwg(1:101))
plot(w(1:101), exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*cos((w(1:101)-0.06)*pi/0.03))
plot(w(1:101), exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.03))
figure
plot(t, exp(-(t-500).^2/(2*150^2)).*cos(0.06*t))
plot(t, exp(-(t-500).^2/(2*100^2)).*cos(0.06*t))
plot(t, exp(-(t-500).^2/(2*100^2)).*cos(0.06*(t-500)))
genv = exp(-(t-500).^2/(2*100^2)).*cos(0.06*(t-500));
genvw=dctI(genv)/dctfactor;
dctfactor = 1e3/(sqrt(5e3*pi))
genvw=dctI(genv)/dctfactor;
genvw=dctI(genv)*dctfactor;
figure
plot(w, genvw)
plot(w(1:1001), genvw(1:1001))
plot(w(1:101), genvw(1:101))
figure
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.03))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.03))/dctfactor)
2*pi/100
2*pi/102
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.06))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi*1e3))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi*1e3/2))/dctfactor)
1/0.075
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
figure
plot(t, fieldtg)
plot(w(1:1001), evawtg(1:1001))
plot(w(1:1001), evawg(1:1001))
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 10*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
plot(t, fieldtg)
[Jg, J1g, J2g, Jorthg, Jpnormg]
sqnorm(psig(:,end))
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 3*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
figure
plot(w(1:1001), evawg(1:1001))
13*0.06
figure
plot(t, evat)
plot(t, evatg)
plot(t, evmiutg)
plot(t, fieldtg)
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(w(1:1001), evaw(1:1001))
figure
plot(t, fieldtg)
plot(t, fieldt)
hold on
plot(t, fieldtg, 'r')
figure
plot(w(1:101), fieldwg(1:101))
hold on
plot(w(1:101), fieldw(1:101), 'r')
whos
[fi0, E0, x, E, P240] = gsV(Vabs240, [-240 240], Nx240);
max(abs(fi0-fi0240))
E-E240
[fi0, E0, x, E, P240] = gsV(Vabs240, [-240 240], Nx240)
x240(1)
x240(end)
max(abs(x-x240))
figure
plot(x240, real(Vabs240))
figure
plot(x240, conj(P240(:, 200)).*P240(:,200))
plot(x240, conj(P240(:, 768)).*P240(:,768))
[fi0, E0, x, E, P240] = gsV(Vabs240, [-240 240], Nx240);
E
[E E240]
whos
figure
plot(x240, V0abs)
plot(x240, V0240)
[fi0, E0, x, E, P240] = gsV(V0240, [-240 240], Nx240);
max(abs(fi0-fi0240))
max(abs(E-E240))
[fi0, E0, x, E, P240] = gsV(Vf, [-240 240], Nx240);
whos
Vf = @(x)1-1./sqrt(x.^2+1)
[fi0, E0, x, E, P240] = gsV(Vf, [-240 240], Nx240);
max(abs(E-E240))
clear fi0 E0 x E
psiE = P240\psi;
figure
psiEg = P240\psig;
plot(t, conj(psiEg).*psiEg)
plot(t, conj(psiE).*psiE)
E(4)-E(1)
E240(4)-E240(1)
E240(2:5)-E0240
1-E0240
[fieldt78g, fieldw78g, psi78g, evat78g, evaw78g, evmiut78g, evmiuw78g, mniterc78g, J78g, J178g, J278g, Jorth78g, Jpnorm78g] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 2e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[J78g, J178g, J278g, Jorth78g, Jpnorm78g]
sqnorm(psig(:,end))
sqnorm(psi(:,end))
[fieldt78g, fieldw78g, psi78g, evat78g, evaw78g, evmiut78g, evmiuw78g, mniterc78g, J78g, J178g, J278g, Jorth78g, Jpnorm78g] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[J78g, J178g, J278g, Jorth78g, Jpnorm78g]
figure
plot(w(1:1001), evaw78g(1:1001))
[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, relE78, conv78, niter78, mallniterc78, J178, maxgrad78, weight78] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldt78g)
hold on
plot(t, fieldt78, 'r')
sqnorm(psi78(:,end))
sqnorm(psi(:,end))
plot(t, fieldt, 'g')
figure
plot(w(1:1001), evaw78g(1:1001))
hold on
plot(w(1:1001), evaw78(1:1001), 'r')
figure
plot(w(1:101), fieldwg(1:101))
hold on
plot(w(1:101), fieldw78(1:101), 'r')
[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, mniterc78, J78, J178, J278, Jorth78, Jpnorm78] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw78, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
plot(w(1:101), fieldw78(1:101), 'r')
plot(w(1:101), fieldw78(1:101), 'g')
plot(w(1:1001), evaw78(1:1001), 'g')
[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, mniterc78, J78, J178, J278, Jorth78, Jpnorm78] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw78, @(w) 1e2*0.5*(1-tanh(100*(w-0.07))), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, relE78, conv78, niter78, mallniterc78, J178, maxgrad78, weight78] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sqnorm(psi78(:,end))
[fieldt78g, fieldw78g, psi78g, evat78g, evaw78g, evmiut78g, evmiuw78g, mniterc78g, J78g, J178g, J278g, Jorth78g, Jpnorm78g] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
figure
plot(t, fieldt78g)
hold on
plot(t, fieldt78, 'r')
figure
plot(w(1:101), fieldwg(1:101))
hold on
plot(w(1:101), fieldw78(1:101),'r')
figure
plot(w(1:1001), evaw78g(1:1001))
hold on
plot(w(1:1001), evaw78(1:1001), 'r')
figure
plot(w(1:1001), evaw78g(1:1001))
figure
plot(w(1:1001), evmiuw78(1:1001))
save field78 fieldt78 fieldw78 evat78 evaw78 evmiut78 evmiuw78 conv78 niter78 J178
whos
%-- 18/08/2014 07:45 --%
load coulomb_optV240
load fieldw78
whos
[~, ~, ~, ~, P240] = gsV(Vf, [-240 240], Nx240);
Vf = @(x)1-1./sqrt(x.^2+1)
[~, ~, ~, ~, P240] = gsV(Vf, [-240 240], Nx240);
fieldiw = instwcos(fieldt78, 1e3);
whos
clear fieldw
load field78
[fieldt78, fieldw78, psi78, evat78, evaw78, evmiut78, evmiuw78, mniterc78, J78, J178, J278, Jorth78, Jpnorm78] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw78, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
conv78(end)-J78
psiE78 = P240\psi78;
figure
plot(t, conj(psiE78).*psiE78)
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(t, conj(psiE78).*psiE78)
E240(1:8)-E0240
figure
plot(t, fieldt78)
fieldiw = instwcos(fieldt78, 1e3);
figure
plot(t, fieldiw)
figure
plot(t, evat78)
figure
plot(t, evmiut78)
npsi78=sqnorm(psi78);
figure
plot(t, npsi78)
evaiw = instw(evat78, 1e3);
figure
plot(t, evaiw)
evmiuiw = instw(evmiut78, 1e3);
figure
plot(t, evmiuiw)
evmiuiw = instw(evmiut78, 0.2);
evaiw = instw(evat78, 0.2);
figure
plot(t, evaiw)
dctfactor = 1e3/(sqrt(5e3*pi))
evawh = dctI(fieldt78(2501:end))*dctfactor/sqrt(2);
figure
plot(w(1:1001), evawwh(1:1001))
plot(w(1:1001), evawh(1:1001))
size(evawh)
wh=0:pi/0.5e3:pi/0.2;
plot(wh(1:1001), evawh(1:1001))
evawh = dctI(evat78(2501:end))*dctfactor/sqrt(2);
plot(wh(1:1001), evawh(1:1001))
evawh1 = dctI(evat78(1:2501))*dctfactor/sqrt(2);
figure
plot(wh(1:1001), evawh1(1:1001))
plot(wh(1:501), evawh1(1:501))
figure
plot(w(1:1001), evaw78(1:1001))
plot(wh(1:501), evawh(1:501))
[fieldt39, fieldw39, psi78, evat39, evaw39, evmiut39, evmiuw39, mniterc39, J39, J139, J239, Jorth39, Jpnorm39] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
[J39, J139, J239, Jpnorm39]
[fieldt39, fieldw39, psi78, evat39, evaw39, evmiut39, evmiuw39, mniterc39, J39, J139, J239, Jorth39, Jpnorm39] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
[J39, J139, J239, Jpnorm39]
figure
plot(w(1:1001), evaw39(1:1001))
%-- 19/08/2014 14:06 --%
load coulomb_optV240
Vf = @(x)1-1./sqrt(x.^2+1)
[~, ~, ~, ~, P240] = gsV(Vf, [-240 240], Nx240);
[fieldt39, fieldw39, psi78, evat39, evaw39, evmiut39, evmiuw39, mniterc39, J39, J139, J239, Jorth39, Jpnorm39] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
[J39, J139, J239, Jpnorm39]
[fieldt39, fieldw39, psi39, evat39, evaw39, evmiut39, evmiuw39, relE39, conv39, niter39, mallniterc39, J139, maxgrad39, weight39] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e-4, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
t=0:0.2:1e3;
wh=0:pi/0.5e3:pi/0.2;
plot(t, fieldt39)
[fieldt39, fieldw39, psi39, evat39, evaw39, evmiut39, evmiuw39, relE39, conv39, niter39, mallniterc39, J139, maxgrad39, weight39] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldt39)
figure
plot(w(1:1001), evaw39(1:1001))
w=0:pi/1e3:pi/0.2;
plot(w(1:1001), evaw39(1:1001))
save field39 fieldt39 fieldw39 evat39 evaw39 evmiut39 evmiuw39 conv39 niter39 J139
figure
plot(w(1:101), fieldw39(1:101))
%-- 24/08/2014 07:52 --%
load coulomb_optV240
whos
Vf = @(x)1-1./sqrt(x.^2+1)
[fieldt39, fieldw39, psi78, evat39, evaw39, evmiut39, evmiuw39, mniterc39, J39, J139, J239, Jorth39, Jpnorm39] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw39, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
whos
load field39_2
[fieldt39, fieldw39, psi78, evat39, evaw39, evmiut39, evmiuw39, mniterc39, J39, J139, J239, Jorth39, Jpnorm39] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.1*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw39, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
[J39, J139, J239, Jpnorm39]
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
plot(t, fieldt39)
figure
plot(w(1:101), fieldw39(1:101))
hold on
plot(w(1:101), exp(-(w(1:101)-0.06).^2/(2*0.01^2)),'r')
figure
plot(w(1:1001), evaw39(1:1001))
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw39, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 100*rectanglefun(w, 0.38, 0.40), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw39, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 100*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
figure
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw39, @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 200*rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldt54)
figure
plot(w(1:101), fieldw54(1:101))
figure
plot(w(1:1001), evaw54(1:1001))
fieldiw = instwcos(fieldt54, 1e3);
figure
plot(t, fieldiw)
save field54 fieldt54 fieldw54 evat54 evaw54 evmiut54 evmiuw54 conv54 niter54 J154
%-- 26/08/2014 12:10 --%
load coulomb_optV240
figure
t=0:0.2:1e3;
w=0:pi/1e3:pi/0.2;
fieldwg = exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015);
plot(w(1:1001), fieldwg)
plot(w(1:1001), fieldwg(1:1001))
plot(w(1:101), fieldwg(1:101))
dctfactor = 1e3/(sqrt(5e3*pi))
fieldtg = dctI(fieldwg)/dctfactor;
plot(t, fieldtg)
fieldwg1 = exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015);
plot(w(1:101), fieldwg1(1:101))
fieldtg1 = dctI(fieldwg1)/dctfactor;
plot(t, fieldtg1)
figure
plot(t, dctI(0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015)/dctfactor)
plot(t, dctI(0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075)))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075))/dctfactor)
plot(t, dctI(exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.00325))/dctfactor)
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 2000*rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 2000*rectanglefun(w, 0.53, 0.55), 0, 1e-4, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 200*rectanglefun(w, 0.53, 0.55), 0, 1e-3, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldt54)
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 200*rectanglefun(w, 0.53, 0.55), 0, 1e-2, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(t, fieldt)
plot(0:0.2:1e3, fieldt)
plot(0:0.2:1e3, allfield(1:6:end))
figure
plot(0:pi/1e3:pi/10, fieldw(1:101))
plot(0:pi/1e3:pi, evaw(1:1001))
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 100*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 2000*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.0075), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 2000*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 3000*rectanglefun(w, 0.53, 0.55), 0, 1e-4, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.1*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 4000*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.5*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 100*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.5*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.53, 0.55), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
figure
plot(w(1:1001), evawg(1:1001))
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.5*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.53, 0.55), 0, 1e-4, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sqnorm(psi54(:,end))
figure
plot(t, fieldt54)
hold on
plot(t, fieldtg, 'r')
figure
plot(w(1:1001), evaw54(1:1001))
Vf = @(x)1-1./sqrt(x.^2+1)
[~, ~, ~, ~, P240] = gsV(Vf, [-240 240], Nx240);
psiE54 = P240\psi54;
figure
%-- 27/08/2014 13:58 --%
prefdir
startdir
pwd
type startup
type startup.m
cd(prefdir)
dir
%-- 27/08/2014 14:06 --%
load coulomb_optV240
[fieldt54, fieldw54, psi54, evat54, evaw54, evmiut54, evmiuw54, relE54, conv54, niter54, mallniterc54, J154, maxgrad54, weight54] = OCfpnorm_evaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.0075) + 0.5*exp(-(w-0.18).^2/(2*0.01^2)).*sin((w-0.18)*pi/0.015), @(w) 1e5*exp(-(w-0.06).^2/(2*0.01^2)) + 1e4*exp(-(w-0.18).^2/(2*0.01^2)), @(w) 10*rectanglefun(w, 0.53, 0.55), 0, 1e-4, 1e3, 0.2, 7, 7, 1e-3, 1e3);
save field54uv fieldt54 fieldw54 evat54 evaw54 evmiut54 evmiuw54 relE54 conv54 niter54 mallniterc54 J154 maxgrad54 weight54
%-- 29/08/2014 08:02 --%
