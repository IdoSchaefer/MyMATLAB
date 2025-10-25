function [Fi, Hts, field] = ihfieldMLS(psi, miu, tsi, H0, allfield)
    [dim, Nt_ts] = size(psi);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    for ti = 1:Nt_ts
        Vtts(:, :, ti) = -miu*allfield((tsi-1)*(Nt_ts - 1) + ti);
    end
    field = allfield((tsi-1)*(Nt_ts - 1) + 1);
    Vhalfv = Vtts(:, :, tmidi);
    Hts = H0 + Vhalfv;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vhalfv)*psi(:, ti);
    end    
end