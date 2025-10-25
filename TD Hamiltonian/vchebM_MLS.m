function M = vchebM_MLS(Hcheb, u0, dim, Ncheb)
    M = zeros(dim, Ncheb);
    M(:, 1) = u0;
    M(:, 2) = Hcheb*u0;
    for k = 3:Ncheb
        M(:, k) = 2*Hcheb*M(:, k-1) - M(:, k-2);
    end
end