dx = 16/128;
x = (-8:dx:7.875).';
fi0 = pi^(-1/4)*exp(-x.^2/2)/sqrt(8);
T = 10; Nt_ts = 6; Ncheb = 15;
Nsamp = 20;
allmv = zeros(1, Nsamp);
aller = zeros(1, Nsamp);
for degi = 1:Nsamp
    deg = log10(200) + (degi-1)*0.1;
    Nt = round(10^deg);
    [U mniter matvecs] = TDHcheb_tsnf(@(x) x.^2/2, @(x,t) x*cos(t), [-10 350], fi0, [-8 8], T, Nt, Nt_ts, Ncheb, 1e-3);
    allmv(degi) = matvecs;
    dt = T/Nt;
    t=0:dt:T;
    mx = evx(U, x);
    if ~isfinite(mx(end))
        display('Error.')
    end
    error = mx - (-0.5*sin(t).*t);
    aller(degi) = max(abs(error));
end
plot(log10(allmv), log10(aller))