function fJforb_penal = get_fJforb_penal(penalforbv)
    fJforb_penal = @(allpsi, Nt_ts, integw) fJforb(allpsi, Nt_ts, integw, penalforbv);
end