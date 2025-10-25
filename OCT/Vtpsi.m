function [Vtts, field] = Vtpsi(psits, x, tsi, allchi, Epenal)
    [Nx, Nt_ts] = size(psits);
    Vtts = zeros(Nx, Nt_ts);
    field = -imag(allchi(:, 1, tsi)'*(x.*psits(:, 1)))/Epenal;
    Vtts(:, 1) = -x.*field;
    for ti = 2:Nt_ts
        Vtts(:, ti) = x.*imag(allchi(:, ti, tsi)'*(x.*psits(:, ti)))/Epenal;
    end
end