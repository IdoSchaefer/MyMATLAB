t = 0:0.2:1e3;
w = 0:pi/1e3:pi/0.2;
dw = pi/1e3;
dctfactor = 1e3/(sqrt(5e3*pi));
load coulomb_optV240
cw = 0.1*sin(0.06*(t-500));
cww = dctI(cw)*dctfactor;
cwwfil = cww.*exp(-(w-0.06).^2/(2*0.01^2));
cwfil = dctI(cwwfil)/dctfactor;
cww_con = fieldw20b(cwwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw).';
cw_con = dctI(cww_con)/dctfactor;
all_evaw = zeros(12, 5001);
allJ1 = zeros(1, 12);
allgsoc = zeros(1, 12);
allsqnormend = zeros(1, 12);
allEend = zeros(1, 12);
for cwi = 1:12
    cww_coni = 0.1*cwi*cww_con;
    [~, ~, psii, ~, all_evaw(cwi, :), ~, ~, ~, ~, allJ1(cwi), ~, ~, ~] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cww_coni, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
    gsend = fi0240'*psii(:, end);
    allgsoc(cwi) = gsend*conj(gsend);
    allsqnormend(cwi) = sqnorm(psii(:, end));
    allEend(cwi) = evH(psii(:, end), @(psi) Hpsi(K240, Vabs240, psi));
end