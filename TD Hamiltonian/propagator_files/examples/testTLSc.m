% The program was written to check the efficiency of Hillel Tal-Ezer's
% propagator, for a TLS.
% You may play with the following parameters:
T = 10; Nt = 100; Nt_ts=5; Ncheb = 3; tol=1e-5;
H0 = [1 0; 0 2];
% The initial state:
ui = [1;0];
% Computing the coupling matrix in all the interior time points in all time
% steps:
tic
Vt = coupM4m(@(t) exp(1i*t), T, Nt, Nt_ts);
toc
tic
[U mniter matvecs] = TDHcheb_tsn(H0, Vt, [0 3], ui, T, Nt, Nt_ts, Ncheb, tol);
toc
% The mean number of iterations, for a time step (should be close to 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which make the most of the
% numerical effort:
matvecs
% Computation of the maximal error - the deviation from the analytical
% result:
dt = T/Nt;
t=0:dt:T;
error = conj(U(2, :)).*U(2, :) - sin(t).^2;
maxer = max(abs(error))