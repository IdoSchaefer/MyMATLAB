function [field, psi, relE, conv, niter, mniterc] = OCcheb1(psi0, target, Vf, Edomain, xdomain, fguess, Epenal, T, dt, Nt_ts, Ncheb, tol)
% fguess is a function handle of the form: @(t).
% Vf is the time independent potential energy. It's a function handle of
% the form: f(x), when x is the row vector of the x grid.
% The space domain is: [-max_x, max_x).
    Nt = T/dt;
    Nx = length(psi0);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    tmidi = round(Nt_ts/2);
    s_ext_i = [1:(tmidi - 1), (tmidi + 1):(Nt_ts + 1)];
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
    %psi = zeros(Nx, Nt+1);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    allmniter = 0;
%tic
    V = Vf(x);
    Vts = zeros(Nx, 1);
    %[allpsi, mniter] = solveOCih2(@ihguess, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-3, fguess, t_ts, dt);
    [~, mniter, ~, ~, psi_history] = SemiGlobalH(@(u, t, v) Hpsi(K, V - x*fguess(t), v),...
        @(u1, t1, u2, t2) x*(fguess(t1) - ones(1, Nt_ts)*fguess(t2)).*u1, 0, [], Edomain, psi0, [0, T], Nts, Nt_ts, Ncheb, tol*1e-3,...
        10, 16, false);
    allmniter = allmniter + mniter;
%toc
    %psiT = allpsi(:, Nt_ts, Nt);
    psiT = psi_history.U(:, allt_lasti);
%     for ti = 1:(Nt+1)
%         lastfield(ti) = fguess((ti-1)*dt);
%     end
    lastfield = fguess(0:dt:T);
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
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        end
%        [allchi, ~, mniter] = solveOCih2(@ihchi2, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-3, allpsi, Epenal);
        [~, mniter, ~, ~, chi_history] = SemiGlobalH(@(u, t, v) Hpsi(K, Vts, v), @Hdiff_op_chi, 0, [], Edomain, chiT, [T, 0], Nts, Nt_ts,...
            Ncheb, tol*1e-3, 10, 16, false);
        allmniter = allmniter + mniter;
%        [allpsi, field, mniter] = solveOC(@Vtpsi, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        %[allpsi, field, mniter] = solveOCih(@ihpsi, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-3, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allpsi, ~, mniter] = solveOCih2(@ihpsi2, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-3, allchi, Epenal);
        allmniter = allmniter + mniter;
        %psiT = allpsi(:, Nt_ts, Nt);
        psiT = allpsi(:, allt_lasti);
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
%     psi(:, 1:Nt) = allpsi(:, 1, :);
%     psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    if nargout>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = target'*psiT;
        conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        conv = conv(:, 1:(niter + 1));
    end
    niter
    mniterc = allmniter/(2*niter + 1)
    
    %%%%%%%% Nested functions %%%%%%%%
    function Hdiff_op_psi(psi, t, psi_mid, tmid)
        tsi = round(t(1)/dt) + 1;
        % There is a problem - there is no test point of chi in the
        % history. We should compute it.
        chi_i = allt_lasti - (tsi - 1)*(Nt_ts - 1):-1:(tsi - 2)*(Nt_ts - 1);
        Vts = chi_psi2field(chi_history.U(chi_i(tmidi), psi_mid)) + V;
        for ti = s_ext_i
            Vti = 
        end
    end

    function newfield = chi_psi2field(chi, psi)
        newfield = -imag(chi'*(x.*psi))/Epenal;
    end
end