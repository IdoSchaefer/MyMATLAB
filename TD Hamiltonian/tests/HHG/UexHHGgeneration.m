load coulomb_optV240 K240 Vabs240 xabs240 fi0240
[Uex_all, ~, matvecs, est_errors] = SemiGlobal1(@(u, t, v) -1i*Hpsi(K240, Vabs240 - xabs240*0.1*sech((t-500)/(170)).^2.*cos(0.06*(t-500)), v), @(u1, t1, u2, t2) 1i*xabs240*0.1*(sech((t1-500)/(170)).^2.*cos(0.06*(t1-500)) - sech((t2-500)/(170))^2*cos(0.06*(t2-500))).*u1, 0, [], [], fi0240, 0:1e3, 3e4, 9, 13, eps);
save UexHHG Uex_all