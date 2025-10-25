load energies3LS
noise = (rand(1,201)-0.5)*2;
fidelity = zeros(1, 5);
fieldw = dctIfrom_ig_sym1(allfield, 10, 0.05*(1-cos((0:5)*pi/6))/2, w(1:6));
t_ts_tp = [t_ts(1:6); (t_ts(4)+t_ts(5))/2];
for orderi = 1:5
    field_noise_tpi = dctIintgrid(fieldw.*(1+noise*10^(-(6-orderi))), 10, t_ts_tp);
    psi_noisei = SemiGlobalHparams(@(psi, field, v) H0*v - field*miu*v, @(psi1, field1, psi2, field2) -miu*(psi1*diag(field1 - field2)), 1, field_noise_tpi, [], EdomainGHz, psi0, [0, 10], 200, 7, 7, 1e-8, 10, 16, [], false);
    fidelity(orderi) = 1 - psi_noisei(2, end).*conj(psi_noisei(2, end));
end
    
