function [Efield, es_field]  = getEfieldH(H0, miu, fieldv, n_levels)

    Nfield = length(fieldv);
    Nhilbert = length(miu);
    Nlevels = length(n_levels);
    Efield = zeros(Nlevels, Nfield);
    es_field = zeros(Nhilbert, Nfield, Nlevels);
    Hfieldi = H0;
    [Pfieldi, D] = eig(Hfieldi);
    Efieldi = diag(D);
    [Efieldi, orderE] = sort(real(Efieldi));
    Pfieldi = Pfieldi(:, orderE);
    Efield(:, 1) = Efieldi(n_levels);
    es_field(:, 1, :) = Pfieldi(:, n_levels);
    es_last_dagger = (Pfieldi(:, n_levels))';
    for fieldi = 2:Nfield
        Hfieldi = H0 - fieldv(fieldi)*miu;
        [Pfieldi, D] = eig(Hfieldi);
        Efieldi = diag(D);
        es_last_Pfieldi = es_last_dagger*Pfieldi;
        es_last_proj = es_last_Pfieldi.*conj(es_last_Pfieldi);
        [~, max_proj_i] = max(es_last_proj, [], 2);
        Efield(:, fieldi) = Efieldi(max_proj_i);
        es_field(:, fieldi, :) = Pfieldi(:, max_proj_i);
        es_last_dagger = (Pfieldi(:, max_proj_i))';
    end
end