function [Vtts, field] = VtchiMLS(chits, miu, tsi, allpsi, Epenal)
    [dim, Nt_ts] = size(chits);
    Vtts = zeros(dim, dim, Nt_ts);
    field = -imag(chits(:, 1)'*miu*allpsi(:, 1, tsi))/Epenal;
    Vtts(:, :, 1) = -miu*field;
    for ti = 2:Nt_ts
        Vtts(:, :, ti) = miu*imag(chits(:, ti)'*miu*allpsi(:, ti, tsi))/Epenal;
    end
end