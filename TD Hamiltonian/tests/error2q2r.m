function [allNt, allmv, aller, max_ers] = error2q2r(T, Nt_ts, Nfm, minNt, Nsamp, Niter, Arnoldi)
% The function computes the data of the error decay with reduction of the
% time-step.
if nargin<6
    Niter = 1;
end
if nargin<7
    Arnoldi = false;
end
allNt = zeros(Nsamp, 1);
allmv = zeros(Nsamp, 1);
aller = zeros(Nsamp, 1);
max_ers.texp = zeros(Nsamp, 1);
max_ers.fU = zeros(Nsamp, 1);
max_ers.conv = zeros(Nsamp, 1);
load qubits2HOtest_problem fieldw2mb2 u02modes Hoperations2modesa
load Uex_2q2rtest Uex
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNt(degi) = Nt;
    allfield_degi = fieldw2allfieldt(fieldw2mb2, T, Nt, Nt_ts);
    if Arnoldi
        [u, ~, matvecs, all_est_er] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield_degi, [], [], u02modes, [0 T], Nt, Nt_ts, Nfm, eps, Niter, 16, [], false);
    else
        [u, ~, matvecs, all_est_er] = SemiGlobalHparams(Hoperations2modesa.psi, Hoperations2modesa.diff_psi, 0, allfield_degi, [], [-15 35], u02modes, [0 T], Nt, Nt_ts, Nfm, eps, Niter, 16, [], false);
%    [u, ~, matvecs, all_est_er] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], Nt, Nt_ts, Nkr, eps, Niter, 20, false);
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