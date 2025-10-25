function mx = evx(U, x)
    szU = size(U);
    Nt = szU(2);
    mx = zeros(1, Nt);
    X = diag(x);
    for ti = 1:Nt
        mx(ti) = exval(X, U(:, ti));
    end
end