function fpsi2chi = get_fchi_gate(targets, limits)
    fpsi2chi = @(psis) chis_gate(psis, targets, limits);
end