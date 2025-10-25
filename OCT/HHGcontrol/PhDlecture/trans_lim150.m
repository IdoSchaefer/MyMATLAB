trans_lim = 0.1*exp(-(t-500).^2./(2*150^2)).*sin(0.06*(t-500));
dctfactor = 1e3/(sqrt(5e3*pi));
trans_limw = dctI(trans_lim)*dctfactor;
dw = pi/1e3;
trans_limw_con = fieldw20b(trans_limw.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw);
trans_lim_con = dctI(trans_limw_con)/dctfactor;