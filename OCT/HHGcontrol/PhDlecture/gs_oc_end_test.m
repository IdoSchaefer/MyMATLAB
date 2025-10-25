t1 = 0:0.1:1e3;
w1 = 0:pi/1e3:pi/0.1;
dw = pi/1e3;
dctfactor1 = 1e3/(sqrt(1e4*pi));
load coulomb_optV240
cw1 = 0.1*sin(0.06*(t1-500));
cww1 = dctI(cw1)*dctfactor1;
cwwfil1 = cww1.*exp(-(w1-0.06).^2/(2*0.01^2));
cwfil1 = dctI(cwwfil1)/dctfactor1;
cww_con1 = fieldw20b(cwwfil1.', 5e5*exp(-(w1.'-0.06).^2/(2*0.01^2)), dw).';
cw_con1 = dctI(cww_con1)/dctfactor1;
all_evaw1 = zeros(12, 10001);
allJ11 = zeros(1, 12);
allgsoc1 = zeros(1, 12);
allsqnormend1 = zeros(1, 12);
allEend1 = zeros(1, 12);
for cwi = 1:12
    cww_coni = 0.1*cwi*cww_con1;
    [~, ~, psii, ~, all_evaw1(cwi, :), ~, ~, ~, ~, allJ11(cwi), ~, ~, ~] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_coni, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.1, 9, 9, 1e-6);
    gsend = fi0240'*psii(:, end);
    allgsoc1(cwi) = gsend*conj(gsend);
    allsqnormend1(cwi) = sqnorm(psii(:, end));
    allEend1(cwi) = evH(psii(:, end), @(psi) Hpsi(K240, Vabs240, psi));
end