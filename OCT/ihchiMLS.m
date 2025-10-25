function [Fi, Hts, field] = ihchiMLS(chits, miu, tsi, H0, allpsi, Epenal)
    [dim, Nt_ts] = size(chits);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    field = -imag(chits(:, 1)'*miu*allpsi(:, 1, tsi))/Epenal;
    Vtts(:, :, 1) = -miu*field;
    for ti = 2:Nt_ts
        Vtts(:, :, ti) = miu*imag(chits(:, ti)'*miu*allpsi(:, ti, tsi))/Epenal;
    end
    Vthalf = Vtts(:, :, tmidi);
    Hts = H0 + Vthalf;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vthalf)*chits(:, ti);
    end    
end