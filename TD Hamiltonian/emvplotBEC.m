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
gs = gsNLHdiag(H, @(u,x,t) conj(u).*u, x, 2e-12);
T = 10; Nt_ts = 5; Ncheb = 12;
minNt = 200;
Niter1st = 4; Niterextrap = 2;
% T = 10; Nt_ts = 7; Ncheb = 7;
% minNt = 350;
Nsamp = 15;
allmv = zeros(1, Nsamp);
aller = zeros(1, Nsamp);
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    [U matvecs] = NLHcheckE(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-10 350], exp(1i*x).*gs, [-8 8], T, Nt, Nt_ts, Ncheb, Niter1st, Niterextrap);
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
plot(log10(allmv), log10(aller), 'o')
xlabel('log(matvecs)')
ylabel('log(error)')