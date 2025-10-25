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
gs = gsNLHdiag(H, @(u,x,t) conj(u).*u, x, 1e-12);
T = 10; Nt_ts = 5; Ncheb = 12; %Nt = 1000;
Nsamp = 15;
allmv = zeros(1, Nsamp);
aller = zeros(1, Nsamp);
for degi = 1:Nsamp
    degtol = -(5 + (degi-1)*0.5);
    tol = 10^degtol;
    degNt = 3 + 0.1*(degi - 1);
    Nt = round(10^degNt);
    [U mniter matvecs] = NLHchebn(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-10 350], exp(1i*x).*gs, [-8 8], T, Nt, Nt_ts, Ncheb, tol);
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
plot(log10(allmv), log10(aller))