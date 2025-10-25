function [Fi, Vts, field] = ihguess(psi, x, tsi, Vvec, fieldf, t_ts, dt)
    [Nx, Nt_ts] = size(psi);
    Vtts = zeros(Nx, Nt_ts);
    Fi = zeros(Nx, Nt_ts);
    tmidi = round(Nt_ts/2);
    for ti = 1:Nt_ts
        Vtts(:, ti) = -x*fieldf((tsi-1)*dt + t_ts(ti));
    end
    field = fieldf((tsi-1)*dt);
    Vhalfv = Vtts(:, tmidi);
    Vts = Vvec + Vhalfv;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, ti) - Vhalfv).*psi(:, ti);
    end    
end