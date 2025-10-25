function M = vchebM(K, V, u0, dim, leftb, rightb, Ncheb)
    M = zeros(dim, Ncheb);
    M(:, 1) = u0;
    M(:, 2) = Hpsicheb(K, V, u0, leftb, rightb);
    for k = 3:Ncheb
        M(:, k) = 2*Hpsicheb(K, V, M(:, k-1), leftb, rightb) - M(:, k-2);
    end
end