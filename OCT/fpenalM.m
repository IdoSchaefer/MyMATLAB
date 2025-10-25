function beta = fpenalM(penalfun, T, dt, Nt_ts)
% penalfun is a function handle of w.
    tcheb = -cos(((1:(Nt_ts - 1)) - 1)*pi/(Nt_ts-1));
    % Note that here, the right boundary is not included - there are
    % Nt_ts-1 points.
    t_ts = 0.5*(tcheb+1)*dt;
    W = pi/dt;
    Nt = T/dt;
    timegrid = zeros(1, Nt*(Nt_ts - 1) + 1);
    for ti = 1:(Nt_ts - 1)
        timegrid(ti:(Nt_ts - 1):(Nt*(Nt_ts - 1) + 1)) = t_ts(ti):dt:T;
    end
    beta = zeros(Nt*(Nt_ts - 1) + 1);
    betafun = @(w, ti, tj) penalfun(w).*cos(w*timegrid(ti)).*cos(w*timegrid(tj));
    for ti = 1:(Nt*(Nt_ts - 1) + 1)
        for tj = 1:(ti - 1)
%            beta(ti, tj) = quadgk(@(w) penalfun(w).*cos(w*timegrid(ti)).*cos(w*timegrid(tj)), 0, W);
            beta(ti, tj) = quadgk(@(w) betafun(w, ti, tj), 0, W);
            beta(tj, ti) = beta(ti, tj);
        end
%        beta(ti, ti) = quadgk(@(w) penalfun(w).*cos(w*timegrid(ti)).*cos(w*timegrid(tj)), 0, W);
        beta(ti, ti) = quadgk(@(w) betafun(w, ti, ti), 0, W);
    end
end
