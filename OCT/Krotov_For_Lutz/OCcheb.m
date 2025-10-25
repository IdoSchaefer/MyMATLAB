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
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    field = zeros(1, Nt+1);
    %lastfield = zeros(1, Nt+1);
    tchebgrid = [kron(ones(1, Nt), t_ts(1:(Nt_ts - 1))) + kron((0:(Nt - 1))*dt, ones(1, Nt_ts - 1)), T];
    integw = chebgridw(Nt, Nt_ts, dt);
    allfield = fguess(tchebgrid);
    lastfield = allfield(1:(Nt_ts - 1):allt_lasti);
%    chi = zeros(Nx, Nt+1);
    psi = zeros(Nx, Nt+1);
    allmniter = 0;
%tic
%    [allpsi, lastfield, mniter] = solveOCg(K, Vf, fguess, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol);
%    [allpsi, lastfield, mniter] = solveOC(@Vtguess, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt);
    [allpsi, ~, mniter] = solveOCih(@ihguess, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-3, fguess, t_ts, dt);
    allmniter = allmniter + mniter;
%toc
    psiT = allpsi(:, Nt_ts, Nt);
%     for ti = 1:(Nt+1)
%         lastfield(ti) = fguess((ti-1)*dt);
%     end
    if nargout>3
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
        if nargout>3
%            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum([0.5*field(1)^2, field(2:Nt).^2, 0.5*field(Nt + 1)^2])*dt;
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(integw.*allfield.^2);
        end
%        [allchi, field1, mniter] = solveOC(@Vtchi, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-2, allpsi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allchi, ~, mniter] = solveOCih(@ihchi, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-3, allpsi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        allmniter = allmniter + mniter;
%        [allpsi, field, mniter] = solveOC(@Vtpsi, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allpsi, field(1:Nt), mniter] = solveOCih(@ihpsi, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-3, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        for tsi = 1:Nt
            allfield((tsi - 1)*(Nt_ts - 1) + 1) = field(tsi);
            for ti = 2:(Nt_ts - 1)
                allfield((tsi - 1)*(Nt_ts - 1) + ti) = -imag(allchi(:, Nt_ts - ti + 1, Nt - tsi + 1)'*(x.*allpsi(:, ti, tsi)))/Epenal;
            end
        end
        allfield(allt_lasti) = -imag(allchi(:, 1, 1)'*(x.*allpsi(:, Nt_ts, Nt)))/Epenal;
        field(Nt + 1) = allfield(allt_lasti);
%         field(1) = -imag(allchi(:, Nt_ts, Nt)'*(x.*allpsi(:, 1, 1)))/Epenal;
%         for ti = 2:Nt
%             field(ti) = -imag(allchi(:, 1, Nt - ti + 2)'*(x.*allpsi(:, 1, ti)))/Epenal;
%         end
        %norm(field(1:Nt)-field1)/norm(field1)
        %field(Nt + 1) = -imag(allchi(:, 1, 1)'*(x.*allpsi(:, Nt_ts, Nt)))/Epenal;
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
    if nargout>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = target'*psiT;
        %conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum([0.5*field(1)^2, field(2:Nt).^2, 0.5*field(Nt + 1)^2])*dt;
        conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(integw.*allfield.^2);
        conv = conv(:, 1:(niter + 1));
    end
    niter
    mniterc = allmniter/(2*niter + 1)
end