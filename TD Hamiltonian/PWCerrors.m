function all_ers = PWCerrors(Uts, Hop, params, Edomain, t_interval, Ncheb)
    t0 = t_interval(1);
    tf = t_interval(2);
    T = tf - t0;
    Nt = size(Uts, 2) - 1;
    dt = T/Nt;
    all_ers = zeros(1, Nt);
    for ti = 1:Nt
        uPWC = SchrPWCcheb(Hop, Uts(:, ti), params(:, ti), Edomain, dt, 1, Ncheb);
        %uPWC = fMchebop(Gop, ev_domain(2), ev_domain(1), @(G) exp(G), Uts(:, ti), dt, 1, Ncheb);
        all_ers(ti) = norm(uPWC(:,2) - Uts(:, ti + 1))/norm(Uts(:, ti + 1));
    end
end