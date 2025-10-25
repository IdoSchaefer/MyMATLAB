function [Vtts, field] = Vtguess(psi, x, tsi, fieldf, t_ts, dt)
    [Nx, Nt_ts] = size(psi);
    Vtts = zeros(Nx, Nt_ts);
    for ti = 1:Nt_ts
        Vtts(:, ti) = -x*fieldf((tsi-1)*dt + t_ts(ti));
    end
    field = fieldf((tsi-1)*dt);
end
