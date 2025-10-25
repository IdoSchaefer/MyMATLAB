load coulomb_optV240
gs_oc_endw_sec3
[~, ~, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
Vf = @(x) 1 - 1./sqrt(1+x.^2);
[~, ~, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
tmidi = round(Nt_ts/2)
Nt_ts=7
tmidi = round(Nt_ts/2)
test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2
gs_oc_endw_sec3
test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2
gs_oc_end_sec3
cw_con_allt_sec_tp = dctIintgrid(cwwsec_con, T, [t_ts(1:(Nt_ts-1)), test_tpoint)/dctfactor;
cw_con_allt_sec_tp = dctIintgrid(cwwsec_con, T, [t_ts(1:(Nt_ts-1)), test_tpoint])/dctfactor;
size(cw_con_allt_sec_tp)
size(cw_con_allt_sec)
[U, mniter, matvecs, max_errors] = SemiGlobalHparams(@(psi, field, v) (H03 - miu240_3*field)*v , @(psi1, field1, psi2, field2) -miu240_3*psi1*(field1.' - field2), 1, cw_con_allt_sec_tp, [], [0 1], [1;0;0], 0:0.2:1e3, 5e3, 7, 7, 1e-5);