% The program demonstrates the application of TDHgeneral, for a TLS.
% You may play with the following parameters:
T = 10; Nt = 100; Nt_ts=5; Ncheb = 3; tol=1e-5;
% The initial state:
ui = [1;0];
H0 = @(v) [1 0;
           0 2]*v;
Vt = @(v, u, t) [0            exp(1i*t);
                 exp(-1i*t)   0         ]*v;
tic
[U mniter matvecs] = TDHgeneral(H0, Vt, [], [0 3], ui, [0 T], Nt, Nt_ts, Ncheb, tol);
toc
% The mean number of iterations, for a time step (should be 1, for ideal
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