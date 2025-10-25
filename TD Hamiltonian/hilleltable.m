% dx = 16/128;
% x = (-8:dx:7.875).';
% % The potential energy matrix:
% V = diag(x.^2/2);
% dp = 2*pi/16;
% p = (-pi/dx):dp:(pi/dx - dp);
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
H = K + V;
gs = gsNLHdiag(H, @(u,x,t) conj(u).*u, x, 2e-12);
T = 10; Nt_ts = 9; Ncheb = 9;
minNt = 140;
Niter1st = 3; Niterextrap = 1;
% T = 10; Nt_ts = 7; Ncheb = 12;
% minNt = 200;
Nsamp = 20;
allNt = zeros(Nsamp, 1);
allmv = zeros(Nsamp, 1);
aller = zeros(Nsamp, 1);
load exactBECu10_2
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
%    deg = log10(minNt) + (degi-1)*0.02;
    Nt = round(10^deg);
    allNt(degi) = Nt;
    [U matvecs] = NLHcheckE(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-1 188], exp(1i*8*x).*gs, [-xdlength/2 xdlength/2], T, Nt, Nt_ts, Ncheb,...
        Niter1st, Niterextrap);
    allmv(degi) = matvecs;
    if ~isfinite(U(1, end))
        display('Error.')
    end
    error = norm(U(:, end) - exactBECu10);
    aller(degi) = error;
end
figure
plot(log10(allmv), log10(aller), 'o')
xlabel('log(matvecs)')
ylabel('log(error)')
data = [allNt allmv aller];