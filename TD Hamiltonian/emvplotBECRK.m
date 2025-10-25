dx = 16/128;
x = (-8:dx:7.875).';
% The potential energy matrix:
V = diag(x.^2/2);
dp = 2*pi/16;
p = ((-pi/dx):dp:(pi/dx - dp)).';
% The kinetic energy matrix in the p domain:
Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
K = 128*ifft(ifft(ifftshift(Kp))')';
% The Hamiltonian:
H = K + V;
gs = gsNLHdiag(H, @(u,x,t) conj(u).*u, x, 2e-12);
p1 = (0:(2*pi/16):(2*pi*(1/dx - 1/16))).';
p1(65:128) = p1(65:128) - 2*pi/dx;
Kv = p1.^2/2;
T = 10;
Nsamp = 15;
allmv = zeros(1, Nsamp);
aller = zeros(1, Nsamp);
for degi = 1:Nsamp
    deg = log10(1500) + (degi-1)*0.1;
    Nt = round(10^deg);
    dt = T/Nt;
    U = RK4(@SEnl, [0 T], exp(1i*x).*gs, dt, Kv, @(u, x, t) x.^2/2 + conj(u).*u, x);
    matvecs = Nt*4;
%    [U mniter matvecs] = NLHchebn(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-10 350], exp(1i*x).*gs, [-8 8], T, Nt, Nt_ts, Ncheb, 1e-3);
    allmv(degi) = matvecs;
    dt = T/Nt;
    t=0:dt:T;
    mx = evx(U, x);
    if ~isfinite(mx(end))
        display('Error.')
    end
    error = mx - sin(t);
    aller(degi) = max(abs(error));
end
figure
plot(log10(allmv), log10(aller), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')