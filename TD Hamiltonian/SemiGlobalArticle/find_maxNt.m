function [umax, maxNt, ermax, allNt, allmv, aller] = find_maxNt(T, Nt_ts, Nkr, minNt)

load coulomb_optV240
uold = TDHxpKr1(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech(-(t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], minNt, Nt_ts, Nkr, eps, 1, 20, false);
deg = log10(minNt) + 0.1;
Nt = round(10^deg);
[u, ~, matvecs] = TDHxpKr1(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech(-(t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], Nt, Nt_ts, Nkr, eps, 1, 20, false);
allmv = matvecs;
allNt = Nt;
old_error = Inf;
error = norm(uold(:, 2) - u(:, 2))/norm(u(:, 2));
aller = error;
degi = 2;
while error<old_error
    old_error = error;
    uold = u;
    oldNt = Nt;
    deg = log10(Nt) + degi*0.1;
    Nt = round(10^deg);
    [u, ~, matvecs] = TDHxpKr1(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech(-(t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], Nt, Nt_ts, Nkr, eps, 1, 20, false);
    if ~isfinite(u(:,2))
        display('Error.')
    end
    error = norm(uold(:, 2) - u(:, 2))/norm(u(:, 2));
    allNt = [allNt; Nt];
    allmv = [allmv; matvecs];
    aller = [aller; error];
    degi = degi + 1;
end
maxNt = oldNt;
umax = uold(:, 2);
ermax = old_error;
figure
plot(log10(allmv), log10(aller), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')

end