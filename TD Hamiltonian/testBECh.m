% The program assumes the existence of the variables:
% T, Nt, Nt_ts, Ncheb, Niter1st, Niterextrap, gs.
dx = sqrt(2*pi/128);
xdlength = 128*dx;
x = (-xdlength/2:dx:(xdlength/2 - dx)).';
%fi0 = pi^(-1/4)*exp(-x.^2/2)/sqrt(8);
load exactBECu10_2
tic
%[U mniter matvecs] = NLHchebn(@(x) x.^2/2, @(u, x, t) conj(u).*u + x*cos(t), [-10 350], gs, [-8 8], T, Nt, Nt_ts, Ncheb, tol);
[U matvecs] = NLHcheckE(@(x) x.^2/2, @(u, x, t) conj(u).*u, [-1 188], exp(1i*8*x).*gs, [-xdlength/2 xdlength/2], T, Nt, Nt_ts, Ncheb, ...
    Niter1st, Niterextrap);
toc
matvecs
if ~isfinite(U(1, end))
    display('Error.')
end
error = norm(U(:, end) - exactBECu10)