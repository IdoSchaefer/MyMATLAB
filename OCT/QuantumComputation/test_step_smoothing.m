load HTLSs_harmonic
Nsigma = 20;
Dsigma = 0.1;
allJmax_overlap4c_smooth = zeros(1, Nsigma);
for sigmai = 1:Nsigma
    sigma = Dsigma*(2*pi/50)*sigmai;
    fEsmoothi = @(t) soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), t, @(x, a, b) softrectfun_erf(x, a, b, sigma));
    psi4c_smoothi = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEsmoothi, [], [-15 15], u0_3, [0, 6], 160, 9, 7,...
        1e-6, 10, 16, [], false);
    allJmax_overlap4c_smooth(sigmai) = Uoverlap_gate(psi4c_smoothi(:, end), target3LSu(:), [1 4 7]);
end
figure
plot((1:Nsigma)*Dsigma, log10(abs(allJmax_overlap4c_smooth - Jmax_overlap4cPWC20)))