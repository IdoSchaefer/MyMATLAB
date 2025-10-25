function [allNt, allmv, aller, max_ers] = error_forced_harmonic(Nt_ts, Nfm, minNt, Nsamp, Niter, Arnoldi)
% The function computes the data of the error decay with reduction of the
% time-step.
if nargin<5
    Niter = 1;
end
if nargin<6
    Arnoldi = false;
end
allNt = zeros(Nsamp, 1);
allmv = zeros(Nsamp, 1);
aller = zeros(Nsamp, 1);
max_ers.texp = zeros(Nsamp, 1);
max_ers.FU = zeros(Nsamp, 1);
max_ers.conv = zeros(Nsamp, 1);
load Uexact_forced_harmonic Uex
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
% The harmonic oscillator ground state:
fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNt(degi) = Nt;
    if Arnoldi
        [u, ~, matvecs, all_est_er] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, [0 10], Nt, Nt_ts, Nfm, eps, Niter, 16, false);
    else
        [u, ~, matvecs, all_est_er] = SemiGlobal(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v), @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, [0 10], Nt, Nt_ts, Nfm, eps, Niter, 16, false);
    end
    allmv(degi) = matvecs;
    max_ers.texp(degi) = all_est_er.texp;
    max_ers.fU(degi) = all_est_er.fU;
    max_ers.conv(degi) = all_est_er.conv;
    if ~isfinite(u(:,2))
        fprintf('\nError.\n')
    end
    error = norm(u(:, 2) - Uex(:, end))/norm(Uex(:, end));
    aller(degi) = error;
end
figure
plot(log10(allmv), log10(aller), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
end