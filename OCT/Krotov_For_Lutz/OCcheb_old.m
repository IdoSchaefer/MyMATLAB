function [field, psi, relE, conv, niter, mniterc] = OCcheb(psi0, target, Vf, Edomain, xdomain, fguess, Epenal, T, dt, Nt_ts, Ncheb, tol)
% fguess is a function handle of the form: @(t).
% Vf is the time independent potential energy. It's a function handle of
% the form: f(x), when x is the row vector of the x grid.
% The space domain is: [-max_x, max_x).
    Nt = T/dt;
    Nx = length(psi0);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
%    dx = 2*max_x/Nx;
%    x = (-max_x:dx:(max_x - dx)).';
%    p = (0:(pi/max_x):(2*pi/dx - pi/max_x)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
%    field = zeros(1, Nt+1);
%    lastfield = zeros(1, Nt+1);
%    chi = zeros(Nx, Nt+1);
    psi = zeros(Nx, Nt+1);
    allmniter = 0;
%tic
%    [allpsi, lastfield, mniter] = solveOCg(K, Vf, fguess, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol);
%    [allpsi, lastfield, mniter] = solveOC(@Vtguess, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt);
    [allpsi, lastfield, mniter] = solveOCih(@ihguess, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt);
    allmniter = allmniter + mniter;
%toc
    psiT = allpsi(:, Nt_ts, Nt);
%     for ti = 1:(Nt+1)
%         lastfield(ti) = fguess((ti-1)*dt);
%     end
    if nargin>3
        field = lastfield;
    end
%    lastfield = field;
    relE = tol + 1;
    conv = zeros(1, 1001);
    maxNiter = 300;
    niter = 0;
    while relE>tol && niter<maxNiter
%        tic
        overlap = target'*psiT;
        chiT = overlap*target;
        if nargin>3
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        end
%        [allchi, field1, mniter] = solveOC(@Vtchi, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-2, allpsi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allchi, field1, mniter] = solveOCih(@ihchi, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-3, allpsi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        allmniter = allmniter + mniter;
%        [allpsi, field, mniter] = solveOC(@Vtpsi, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allpsi, field, mniter] = solveOCih(@ihpsi, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-3, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        allmniter = allmniter + mniter;
        psiT = allpsi(:, Nt_ts, Nt);
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    if nargin>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = target'*psiT;
        conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        conv = conv(:, 1:(niter + 1));
    end
    niter
    mniterc = allmniter/(2*niter + 1)
end