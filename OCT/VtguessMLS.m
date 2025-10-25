function [Vtts, field] = VtguessMLS(psi, miu, tsi, fieldf, t_ts, dt)
    [dim, Nt_ts] = size(psi);
    Vtts = zeros(dim, dim, Nt_ts);
    for ti = 1:Nt_ts
        Vtts(:, :, ti) = -miu*fieldf((tsi-1)*dt + t_ts(ti));
    end
    field = fieldf((tsi-1)*dt);
end
