T = 1e3;
Nsamp = 17;
minNt = 4e3;
allNt5 = zeros(Nsamp, 1);
allmv5 = zeros(Nsamp, 1);
aller5 = zeros(Nsamp, 1);
all_texp_er5 = zeros(Nsamp, 1);
allFUer5 = zeros(Nsamp, 1);
all_conv_er5 = zeros(Nsamp, 1);
load coulomb_optV240
load Uexact
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNt5(degi) = Nt;
    dt = T/Nt
    [u5, ~, matvecs5, all_est_er5] = TDHxpKr1(K240, Vabs240, @(u,x,t) -xabs240*0.1*sech(-(t-500)/(170)).^2.*cos(0.06*(t-500)), [], fi0240, x240, [0 T], Nt, 5, 5, eps, 1, 3, false);
    allmv5(degi) = matvecs5;
    all_texp_er5(degi) = all_est_er5.texp;
    allFUer5(degi) = all_est_er5.FU;
    all_conv_er5(degi) = all_est_er5.conv;
    if ~isfinite(u5(:,2))
        display('Error.')
    end
    error = norm(u5(:, 2) - Uex(:, end))/norm(Uex(:, end));
    aller5(degi) = error;
end
figure
plot(log10(allmv5), log10(aller5), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
%u5data = [allNt allmv aller];