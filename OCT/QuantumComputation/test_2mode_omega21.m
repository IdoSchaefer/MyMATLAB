load HTLSs_harmonic
N2mode = 101;
all_infidelities2mode_a = zeros(1, N2mode);
all_infidelities2mode_b = zeros(1, N2mode);
H2modesu_i = H2modesu;
isingle2ndmode = [5, 9, 16:18];
is_isingle = false(19, 1);
is_isingle(isingle2ndmode) = true;
is_isingle_diag = diag(is_isingle);
for i2mode = 1:N2mode
    omega21 = i2mode - 1;
    H2modesu_i(is_isingle_diag) = omega21;
    H2modesu_i(19, 19) = 2*omega21;
    Hoperations2modei = Hmats2Hops2(H2modesu_i, H2modes1cu, H2modes2cu);
    [~, ~, psi2modei, ~, ~, ~, ~, Jterms2modei, ~, ~, ~] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modei, 2,...
        fcouplingOp2modes, [-100 170], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3a, 6, 0.025, 9, 9, 1e-6);
    all_infidelities2mode_a(i2mode) = 1 - Jterms2modei.Jmax;
    all_infidelities2mode_b(i2mode) = 1 - Uoverlap_gate(psi2modei(:, end), target2modes, [1 5 9]);
end
figure
plot(0:(N2mode - 1), log10(all_infidelities2mode_a))
hold on
plot(0:(N2mode - 1), log10(all_infidelities2mode_b))
    