t = 0:0.2:1e3;
w = 0:pi/1e3:pi/0.2;
dw = pi/1e3;
dctfactor = 1e3/(sqrt(5e3*pi));
load coulomb_optV240
cwsec = 0.1*sin(0.12*(t-500));
cwwsec = dctI(cwsec)*dctfactor;
cwwfilsec = cwwsec.*exp(-(w-0.12).^2/(2*0.01^2));
cwfilsec = dctI(cwwfilsec)/dctfactor;
cwwsec_con = fieldw20b(cwwfilsec.', 5e5*exp(-(w.'-0.12).^2/(2*0.01^2)), dw).';
cw_consec = dctI(cwwsec_con)/dctfactor;
cw_consec_window = cw_consec*0.5.*(tanh((t-120)/30) - tanh((t-880)/30));
cww_consec_window = dctI(cw_consec_window)*dctfactor;
all_evawsecw = zeros(12, 5001);
allJ1secw = zeros(1, 12);
allgsocsecw = zeros(1, 12);
allsqnormendsecw = zeros(1, 12);
allEendsecw = zeros(1, 12);
for cwi = 1:12
    cwwsec_conwi = 0.1*cwi*cww_consec_window;
    [~, ~, psii, ~, all_evawsecw(cwi, :), ~, ~, ~, ~, allJ1secw(cwi), ~, ~, ~] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), cwwsec_conwi, @(w) 5e5*exp(-(w-0.12).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
    gsendsecw = fi0240'*psii(:, end);
    allgsocsecw(cwi) = gsendsecw*conj(gsendsecw);
    allsqnormendsecw(cwi) = sqnorm(psii(:, end));
    allEendsecw(cwi) = evH(psii(:, end), @(psi) Hpsi(K240, Vabs240, psi));
end