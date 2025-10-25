function [allNt, allmv, aller, allCPU] =  error_decayPWC(Hop, ev_domain, t0tf, u0, uf_exact, Ncheb, minNt, Nsamp)
    allNt = zeros(Nsamp, 1);
    allmv = zeros(Nsamp, 1);
    aller = zeros(Nsamp, 1);
    allCPU = zeros(Nsamp, 1);
    t0 = t0tf(1);
    tf = t0tf(2);
    for degi = 1:Nsamp
        deg = log10(minNt) + (degi-1)*0.1;
        Nt = round(10^deg);
        allNt(degi) = Nt;
        dt = tf/Nt;
        tgrid = ((t0 + dt/2):dt:(tf - dt/2));
        tic
        U = SchrPWCcheb(Hop, u0, tgrid, ev_domain, tf - t0, Nt, Ncheb);
        allCPU(degi) = toc;
        matvecs = Nt*(Ncheb - 1);
        allmv(degi) = matvecs;
        if ~isfinite(U(:, end))
            fprintf('\nError.\n')
        end
        error = norm(U(:, end) - uf_exact)/norm(uf_exact);
        aller(degi) = error;
    end
    figure
    plot(log10(allmv), log10(aller), '-o')
    xlabel('log(matvecs)')
    ylabel('log(error)')
end