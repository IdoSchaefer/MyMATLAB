% The program assumes the existence of the variables:
% T, Nt, Nt_ts, Ncheb, tol, Niter1st, Niterextrap, gs.
dx = 16/128;
x = (-8:dx:7.875).';
%fi0 = pi^(-1/4)*exp(-x.^2/2)/sqrt(8);
tic
%[U mniter matvecs] = NLHchebn(@(x) x.^2/2, @(u, x, t) conj(u).*u + x*cos(t), [-10 350], gs, [-8 8], T, Nt, Nt_ts, Ncheb, tol);
[U matvecs] = NLHcheckE(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-10 350], exp(1i*x).*gs, [-8 8], T, Nt, Nt_ts, Ncheb, Niter1st, Niterextrap);
toc
matvecs
dt = T/Nt;
t=0:dt:T;
mx = evx(U, x);
if ~isfinite(mx(end))
    display('Error.')
end
error = (mx - sin(t));
%error = mx - (-0.5*sin(t).*t);
max(abs(error))