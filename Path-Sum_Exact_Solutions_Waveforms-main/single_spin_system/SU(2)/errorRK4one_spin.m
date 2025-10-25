function [allNt, all_elapsed_times, all_max_errors] = errorRK4one_spin(para, TMAX, minNt, Nsamp)
allNt = zeros(Nsamp, 1);
all_elapsed_times = zeros(Nsamp, 1);
all_max_errors = zeros(Nsamp, 1);
% The first call consumes an extra time, so we do it here:
[~, ~] = one_spin_RK4(para, TMAX, minNt);
for degi = 1:Nsamp
    deg = log10(minNt) + (degi-1)*0.1;
    Nt = round(10^deg);
    allNt(degi) = Nt;
    tic
    for testi = 1:200
        [rho, ~] = one_spin_RK4(para, TMAX, Nt);
    end
    all_elapsed_times(degi) = toc/200;
    rho_exact = one_spin_ODE(para, TMAX, Nt + 1, 1e-13);
    % It is stupid to do it again and again instead of interpolating, but
    % it's a TLS, so who cares?
    f = plotfrob(rho, rho_exact);
    all_max_errors(degi) = max(f);
end
figure
plot(log10(all_elapsed_times), log10(all_max_errors), '-o')
xlabel('log(Number of seconds)')
ylabel('log(error)')
end