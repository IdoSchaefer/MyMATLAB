t1 = 0:0.1:1e3;
w1 = 0:pi/1e3:pi/0.1;
dw = pi/1e3;
dctfactor1 = 1e3/(sqrt(1e4*pi));
load coulomb_optV240
cwsec1 = 0.1*sin(0.12*(t1-500));
cwwsec1 = dctI(cwsec1)*dctfactor1;
cwwfilsec1 = cwwsec1.*exp(-(w1-0.12).^2/(2*0.01^2));
cwfilsec1 = dctI(cwwfilsec1)/dctfactor1;
cwwsec_con1 = fieldw20b(cwwfilsec1.', 5e5*exp(-(w1.'-0.12).^2/(2*0.01^2)), dw).';
cw_consec1 = dctI(cwwsec_con1)/dctfactor1;
all_evawsec1 = zeros(12, 10001);
allJ1sec1 = zeros(1, 12);
allgsocsec1 = zeros(1, 12);
allsqnormendsec1 = zeros(1, 12);
allEendsec1 = zeros(1, 12);
for cwi = 1:12
    cwwsec_coni = 0.1*cwi*cwwsec_con1;
    [~, ~, psii, ~, all_evawsec1(cwi, :), ~, ~, ~, ~, allJ1sec1(cwi), ~, ~, ~] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cwwsec_coni, @(w) 5e5*exp(-(w-0.12).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.1, 9, 9, 1e-6);
    gsend = fi0240'*psii(:, end);
    allgsocsec1(cwi) = gsend*conj(gsend);
    allsqnormendsec1(cwi) = sqnorm(psii(:, end));
    allEendsec1(cwi) = evH(psii(:, end), @(psi) Hpsi(K240, Vabs240, psi));
end