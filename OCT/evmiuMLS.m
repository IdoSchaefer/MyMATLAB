function mmiu = evmiuMLS(psi, miu)
    Nt = size(psi, 2);
    mmiu = zeros(1, Nt);
    for ti = 1:Nt
        mmiu(ti) = psi(:, ti)'*miu*psi(:, ti);
    end
end
