% The program demonstrates the application of TDHgeneral, for a TLS.
% You may play with the following parameters:
T = 100; Nts = 1000; Nt_ts=9; Ncheb = 5; tol=1e-5;
% The output time-grid:
dt = 0.05;
t=0:dt:T;
% The initial state:
ui = [1;
      0];
H0 = [1 0;
      0 2];
mu = [0 1;
      1 0];
Hop = @(u, t, v) H0*v - cos(t)*mu*v;
Hdiff_op = @(u1, t1, u2, t2) -(cos(t1) - cos(t2)).*(mu*u1);   
fprintf('\nSemi-global algorithm computation:\n')
tic
[U, mniter, matvecs] = SemiGlobalH(Hop, Hdiff_op, 1, [], [0, 3], ui, t, Nts, Nt_ts, Ncheb, tol);
toc
% The mean number of iterations, for a time step (should be 1, for ideal
% efficiency):
mniter
% The number of matrix-vector multiplications, which make the most of the
% numerical effort:
matvecs
% Computation of the maximal error:
RKfun = @(t, u) -1i*(H0*u - cos(t)*mu*u);
options = odeset('RelTol', 1e-13, 'absTol', 1e-13);
fprintf('\node45 computation:\n')
tic
[time, URK] = ode45(RKfun, [0 T/2 T], ui, options);
toc
URK = URK(end, :).';
error = norm(U(:, end) - URK(:, end))/norm(URK(:, end))