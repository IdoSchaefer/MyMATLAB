function [Fi, Hts, field] = ihpsiMLS(psits, miu, tsi, H0, allchi, Epenal)
    [dim, Nt_ts] = size(psits);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    field = -imag(allchi(:, 1, tsi)'*miu*psits(:, 1))/Epenal;
    Vtts(:, :, 1) = -miu*field;
    for ti = 2:Nt_ts
        Vtts(:, :, ti) = miu*imag(allchi(:, ti, tsi)'*miu*psits(:, ti))/Epenal;
    end
    Vthalf = Vtts(:, :, tmidi);
    Hts = H0 + Vthalf;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vthalf)*psits(:, ti);
    end    
end