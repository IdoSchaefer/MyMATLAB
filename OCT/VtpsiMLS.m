function [Vtts, field] = VtpsiMLS(psits, miu, tsi, allchi, Epenal)
    [dim, Nt_ts] = size(psits);
    Vtts = zeros(dim, dim, Nt_ts);
    field = -imag(allchi(:, 1, tsi)'*miu*psits(:, 1))/Epenal;
    Vtts(:, :, 1) = -miu*field;
    for ti = 2:Nt_ts
        Vtts(:, :, ti) = miu*imag(allchi(:, ti, tsi)'*miu*psits(:, ti))/Epenal;
    end
end