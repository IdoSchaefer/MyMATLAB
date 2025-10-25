H0 = [1 0;
      0 2];
ui = [1;
      0];
Vt = coupM4(@(t) exp(1i*t), 10, 100, 10);
U = TDHdiag_tsc(H0, Vt, ui, 10, 100, 10, 1e-8);
% The probability:
P = conj(U).*U;
t = 0:0.1:10;
figure
plot(t, P)
% The deviation from the analytic solution:
error = sin(t).^2 - P(2, :);
figure
plot(t, error)
U1 = TDHdiag_tsc(H0, Vt, ui, 10, 100, 10, 1e-9);
U1(:, 1:5)
Ncheb = findNcheb1(@(x) exp(-1i*0.1*x), 1, 2, 1e-10, 2)
U2 = TDHcheb_tsc(H0, 1, 2, Vt, ui, 10, 100, 10, Ncheb, 1e-5);
figure
P2 = conj(U2).*U2;
plot(t, P2)
error2 = sin(t).^2 - P2(2, :);
figure
plot(t, error2)
U3 = TDHcheb_tsc(H0, 1, 2, Vt, ui, 10, 100, 10, Ncheb, 1e-7);
U3(:, 1:5)