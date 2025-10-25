function [allNtRK, allmvRK, allerRK] =  error_forced_harmonicRK4(minNt, Nsamp)
allNtRK = zeros(Nsamp, 1);
allmvRK = zeros(Nsamp, 1);
allerRK = zeros(Nsamp, 1);
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
    allNtRK(degi) = Nt;
    dt = 10/Nt;
    uRK = RK4uf(@(t, u) -1i*Hpsi(K, V + x*cos(t), u), [0 10], fi0, dt);
    matvecs = Nt*4;
    allmvRK(degi) = matvecs;
    if ~isfinite(uRK)
        fprintf('\nError.\n')
    end
    error = norm(uRK - Uex(:, end))/norm(Uex(:, end));
    allerRK(degi) = error;
end
figure
plot(log10(allmvRK), log10(allerRK), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
end