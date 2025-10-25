function [Vtts, field] = Vtchi(chits, x, tsi, allpsi, Epenal)
    [Nx, Nt_ts] = size(chits);
    Vtts = zeros(Nx, Nt_ts);
    field = -imag(chits(:, 1)'*(x.*allpsi(:, 1, tsi)))/Epenal;
    Vtts(:, 1) = -x.*field;
    for ti = 2:Nt_ts
        Vtts(:, ti) = x.*imag(chits(:, ti)'*(x.*allpsi(:, ti, tsi)))/Epenal;
    end
end