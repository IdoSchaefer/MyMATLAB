function [Fi, Vts, field] = ihpsi(psits, x, tsi, Vvec, allchi, Epenal)
    [Nx, Nt_ts] = size(psits);
    Vtts = zeros(Nx, Nt_ts);
    Fi = zeros(Nx, Nt_ts);
    tmidi = round(Nt_ts/2);
    field = -imag(allchi(:, 1, tsi)'*(x.*psits(:, 1)))/Epenal;
    Vtts(:, 1) = -x.*field;
    for ti = 2:Nt_ts
        Vtts(:, ti) = x.*imag(allchi(:, ti, tsi)'*(x.*psits(:, ti)))/Epenal;
    end
    Vhalfv = Vtts(:, tmidi);
    Vts = Vvec + Vhalfv;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, ti) - Vhalfv).*psits(:, ti);
    end    
end