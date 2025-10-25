T = 10; Nt = 100; Npcycle = 5; maxNp = 100; tol = 1e-8;
H = [1 1;
     1 1];
u0 = [1;0];
Hoperation = @(v) H*v;
f = @(x) exp(-1i*x);
[U, Npf, maxEestimated] = fMtLeja(Hoperation, [-1 3], f, u0, T, Nt, Npcycle, maxNp, tol);
dt = T/Nt;
t=0:dt:T;
figure
plot(t, conj(U).*U)
xlabel('$t$ [a.u]', 'interpreter', 'latex')
ylabel('occupation numbers', 'interpreter', 'latex')
error_analytic = conj(U(2, :)).*U(2, :) - sin(t).^2;
maxEanalytic = max(abs(error_analytic))