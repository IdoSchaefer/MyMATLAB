function [mniter, matvecs, psi] = testMLS(psi0, E0, Edomain, miu, fguess, T, dt, Nt_ts, Ncheb, tol)
% The program checks the efficiency of the solver in solving a problem of 3
% level system, using given parameters dt, Nt_ts, Ncheb.
    Nt = T/dt;
    H0 = diag(E0);
    Npsi = length(psi0);
    psi = zeros(Npsi, Nt + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allfield = zeros(1, Nt*(Nt_ts - 1) + 1);
    for tsi = 1:Nt
        for ti = 1:(Nt_ts - 1)
            allfield((tsi - 1)*(Nt_ts - 1) + ti) = fguess((tsi - 1)*dt + t_ts(ti));
        end
    end
    allfield( Nt*(Nt_ts - 1) + 1) = fguess(T);
tic
    [allpsi, field, mniter, matvecs] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol, allfield);
toc
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    if ~isfinite(psi(1, end))
        display('Error')
    end
    mniter
    matvecs
end
