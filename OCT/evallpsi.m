function allev = evallpsi(allpsi, M)
% The function computes the expectation values of M in all the interior
% points within all time steps of allpsi.
% allev is a row vector.
    [Npsi, Nt_ts, Nt] = size(allpsi);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    allev = zeros(1, allt_lasti);
    for tsi = 1:Nt
        for ti = 1:(Nt_ts - 1)
            allev((tsi - 1)*(Nt_ts - 1) + ti) = allpsi(:, ti, tsi)'*M*allpsi(:, ti, tsi);
        end
    end
    allev(allt_lasti) = allpsi(:, Nt_ts, Nt)'*M*allpsi(:, Nt_ts, Nt);
end