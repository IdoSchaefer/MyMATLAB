function [Fi, Vts, field] = ihfield(psi, x, tsi, Vvec, allfield)
    [Nx, Nt_ts] = size(psi);
%    Vtts = zeros(Nx, Nt_ts);
    Fi = zeros(Nx, Nt_ts);
    tmidi = round(Nt_ts/2);
%     for ti = 1:Nt_ts
%         Vtts(:, ti) = -x*allfield((tsi-1)*(Nt_ts - 1) + ti);
%     end
    Vtts = -x*allfield((tsi-1)*(Nt_ts - 1) + (1:Nt_ts));
    field = allfield((tsi-1)*(Nt_ts - 1) + 1);
    Vhalfv = Vtts(:, tmidi);
    Vts = Vvec + Vhalfv;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, ti) - Vhalfv).*psi(:, ti);
    end    
end