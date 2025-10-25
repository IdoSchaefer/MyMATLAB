% The program assumes the existence of the variables:
% T, Nt, Nt_ts, tol.
H0 = [1 0; 0 2];
ui = [1;0];
Vt = coupM4m(@(t) exp(1i*t), T, Nt, Nt_ts);
%Vt = coupM4(@(t) exp(1i*t), T, Nt, Nt_ts);
tic
%[U mniter] = TDHdiag_tsn1(H0, Vt, ui, T, Nt, Nt_ts, tol);
[U mniter] = TDHdiag_tsn2(H0, Vt, ui, T, Nt, Nt_ts, tol);
%U = TDHdiag_tsn1(H0, Vt, ui, T, Nt, Nt_ts, tol);
toc
mniter
%matvecs = mniter*Nt*2*Nt_ts + Nt + Nt_ts + 2
matvecs = mniter*Nt*(3*Nt_ts + 2) + Nt*Nt_ts + Nt_ts + 2
dt = T/Nt;
t=0:dt:T;
error = conj(U(2, :)).*U(2, :) - sin(t).^2;
max(abs(error))