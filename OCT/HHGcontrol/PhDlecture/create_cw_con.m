cw = 0.1*sin(0.06*(t-500));
cww = dctI(cw)*dctfactor;
cwwfil = cww.*exp(-(w-0.06).^2/(2*0.01^2));
cwfil = dctI(cwwfil)/dctfactor;
cww_con = fieldw20b(cwwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw).';
cw_con = dctI(cww_con)/dctfactor;