function [allNtRK, allmvRK, allerRK] =  error2q2rRK(T, minNt, Nsamp)
allNtRK = zeros(Nsamp, 1);
allmvRK = zeros(Nsamp, 1);
allerRK = zeros(Nsamp, 1);
load qubits2HOtest_problem fieldw2mb2 u02modes Hu2modesa H1u2modesa H2u2modesa
load Uex_2q2rtest Uex
dctfactor = T/sqrt(280*pi);
dw = pi/T;
wgrid = 0:dw:215*dw;
fieldw2mb2T = fieldw2mb2(:,1:216).';
get_fieldt = @(t) inv_dctIMtp_wgrid(fieldw2mb2T, t, wgrid, 280).'/dctfactor;
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNtRK(degi) = Nt;
    dt = T/Nt;
    uRK = RK4uf(@(t, u) -1i*(Hu2modesa*u - [H1u2modesa*u, H2u2modesa*u]*get_fieldt(t)), [0 T], u02modes, dt);
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