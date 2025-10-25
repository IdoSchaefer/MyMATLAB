function [Ftaylor, Fdirect, error] = Ftest(H, Tts, Nt_ts)
    [dim, ~]  = size(H);
    [~, D] = eig(H);
    eigval = diag(D);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*Tts;
    Ftaylor = zeros(dim, Nt_ts);
    expeigval = exp(-1i*eigval*t_ts);
    Mtaylor = -1i*eigval*t_ts;
    for ti = 1:Nt_ts
        Mterm = ones(dim, 1);
%        lastMterm = Mterm;
        Ftaylor(:, ti) = Mterm;
        polydeg = 1;
        while norm(Mterm)/norm(Ftaylor(:, ti)) > 1e-16
%        while norm(Mterm)/norm(lastMterm) > 1e-3
%            lastMterm = Mterm;
            Mterm = Mtaylor(:, ti).*Mterm/(polydeg + Nt_ts);
            Ftaylor(:, ti) = Ftaylor(:, ti) + Mterm;
            polydeg = polydeg + 1;
        end
    end
    Ftaylor = Ftaylor.*(ones(dim, 1)*(t_ts.^Nt_ts));
    Fdirect = expeigval;
    Mezer = -1i*eigval*t_ts;
    for polyi = 1:Nt_ts
        Fdirect = polyi*(Fdirect - ones(dim, Nt_ts))./Mezer;
    end
    Fdirect = Fdirect.*(ones(dim, 1)*(t_ts.^Nt_ts));
    error = (Fdirect - Ftaylor)./Ftaylor;
end