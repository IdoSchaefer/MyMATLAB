function [Fi, Hts, field] = ihchigMLS(chi, miu, tsi, H0, fieldf, t_ts, dt, ihterm)
    [dim, Nt_ts] = size(chi);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    for ti = 1:Nt_ts
        Vtts(:, :, ti) = -miu*fieldf((tsi-1)*dt + t_ts(ti));
    end
    field = fieldf((tsi-1)*dt);
    Vhalfv = Vtts(:, :, tmidi);
    Hts = H0 + Vhalfv;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vhalfv)*chi(:, ti) + ihterm(:, ti, tsi);
    end    
end