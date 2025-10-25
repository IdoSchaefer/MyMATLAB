function nU = normUw(U, w)
    szU = size(U);
    Nt = szU(2);
    nU = zeros(1, Nt);
    for ti = 1:Nt
        nU(ti) = sqrt(U(:, ti)'*(w.*U(:, ti)));
    end
end