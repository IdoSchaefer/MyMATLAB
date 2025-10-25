function [allNt, allmv, aller] =  error2q2rPWC(T, Ncheb, minNt, Nsamp)
allNt = zeros(Nsamp, 1);
allmv = zeros(Nsamp, 1);
aller = zeros(Nsamp, 1);
load qubits2HOtest_problem fieldw2mb2 u02modes Hu2modesa H1u2modesa H2u2modesa
load Uex_2q2rtest Uex
dctfactor = T/sqrt(280*pi);
dw = pi/T;
wgrid = 0:dw:215*dw;
fieldw2mb2T = fieldw2mb2(:,1:216).';
Hop = @(v, fields) Hu2modesa*v - [H1u2modesa*v, H2u2modesa*v]*fields;
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNt(degi) = Nt;
    dt = T/Nt;
    tgrid = (dt/2:dt:(T - dt/2)).';
    fieldt = inv_dctIMtp_wgrid(fieldw2mb2T, tgrid, wgrid, 280).'/dctfactor;
    U = SchrPWCcheb(Hop, u02modes, fieldt, [-15 35], T, Nt, Ncheb);
    matvecs = Nt*(Ncheb - 1);
    allmv(degi) = matvecs;
    if ~isfinite(U(:, end))
        fprintf('\nError.\n')
    end
    error = norm(U(:, end) - Uex(:, end))/norm(Uex(:, end));
    aller(degi) = error;
end
figure
plot(log10(allmv), log10(aller), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
end