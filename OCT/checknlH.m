function [psi, evmiut, evmiuw, matvecs, mniter] = checknlH(Vf, psi0, Edomain, xdomain, fieldf, T, dt, Nt_ts, Ncheb, tol)
    Nt = T/dt;
    Nx = length(psi0);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
    psi = zeros(Nx, Nt + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allfield = zeros(1, Nt*(Nt_ts - 1) + 1);
    for tsi = 1:Nt
        for ti = 1:(Nt_ts - 1)
            allfield((tsi - 1)*(Nt_ts - 1) + ti) = fieldf((tsi - 1)*dt + t_ts(ti));
        end
    end
    allfield( Nt*(Nt_ts - 1) + 1) = fieldf(T);
tic
    [allpsi, field, mniter, matvecs] = solveOCih(@ihfield, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol, allfield);
toc
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    if ~isfinite(psi(1, end))
        display('Error')
    end
    evmiut = evmiu(psi, x);
    evmiuw = dctI(evmiut);
    mniter
    matvecs
end    