function [Fi, Vts, field] = ihchi(chits, x, tsi, Vvec, allpsi, Epenal)
    [Nx, Nt_ts] = size(chits);
    Vtts = zeros(Nx, Nt_ts);
    Fi = zeros(Nx, Nt_ts);
    tmidi = round(Nt_ts/2);
    field = -imag(chits(:, 1)'*(x.*allpsi(:, 1, tsi)))/Epenal;
    Vtts(:, 1) = -x.*field;
    for ti = 2:Nt_ts
        Vtts(:, ti) = x.*imag(chits(:, ti)'*(x.*allpsi(:, ti, tsi)))/Epenal;
    end
    Vhalfv = Vtts(:, tmidi);
    Vts = Vvec + Vhalfv;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, ti) - Vhalfv).*chits(:, ti);
    end        
end