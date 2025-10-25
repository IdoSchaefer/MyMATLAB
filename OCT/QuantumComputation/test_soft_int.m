load HTLSs_harmonic
Nsigma = 50;
Dsigma = 0.1;
allJmax_overlap4cInt = zeros(1, Nsigma);
for sigmai = 1:Nsigma
    sigma = Dsigma*(2*pi/50)*sigmai;
    fEInti = @(t) softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), t, @(x, a, b) softrectfun_erf(x, a, b, sigma));
    psi4cInti = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEInti, [], [-15 15], u0_3, [0, 6], 160, 9, 7,...
        1e-6, 10, 16, [], false);
    allJmax_overlap4cInt(sigmai) = Uoverlap_gate(psi4cInti(:, end), target3LSu(:), [1 4 7]);
end
figure
plot((1:Nsigma)*Dsigma, log10(1 - allJmax_overlap4cInt))