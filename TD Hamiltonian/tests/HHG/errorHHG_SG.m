function [allNt, allmv, aller, all_est_ers] = errorHHG_SG(Nt_ts, Nkr, minNt, Nsamp, Niter)
% The function generates the error decay curve for the HHG problem from the
% 2017 semi-global paper.
    if nargin<5
        Niter = 1;
    end
    load coulomb_optV240 K240 Vabs240 xabs240 fi0240
    load Uex_article Uex
    [allNt, allmv, aller, all_est_ers] = error_decaySG2(@(u, t, v) -1i*Hpsi(K240, Vabs240 - xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), v),...
        @(u1, t1, u2, t2) 1i*xabs240*0.1*(sech((t1-500)/(170)).^2.*cos(0.06*(t1-500)) - sech((t2-500)/(170))^2*cos(0.06*(t2-500))).*u1,...
        0, [], [], fi0240, Uex(:, 2), [0 1e3], Nt_ts, Nkr, minNt, Nsamp, Niter);
end