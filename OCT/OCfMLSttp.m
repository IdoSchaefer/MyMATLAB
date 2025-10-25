function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mniterc, J1, maxgrad] = OCfMLSttp(psi0, E0, Edomain, ...
    miu, fguess, penalM, filtermiu, T, dt, Nt_ts, Ncheb, tol)
% fguess is a function handle of the form: @(w).
% E0 is a vector containing the eigenenergies of the time independent
% Hamiltonian matrix.
    Nt = T/dt;
    dw = pi/T;
    dim = length(psi0);
    H0 = diag(E0);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    tolprop = 1e-2*tol;
    fieldw = getwvector(fguess);
    vfiltermiu = getwvector(filtermiu);
    allfield = dctIintgrid(fieldw, T, t_ts(1:(Nt_ts-1))).';
    fieldt = zeros(1, Nt+1);
%    lastfield = zeros(1, Nt+1);
%    chi = zeros(dim, Nt+1);
    psi = zeros(dim, Nt+1);
    chiT = zeros(dim, 1);
    evmiut = zeros(1, Nt + 1);
    evmiuw = zeros(1, Nt + 1);
    evmiufil = zeros(1, Nt + 1);
%     % The small penalty matrix, for computation of the convergence without too much numerical effort, instead of using intchebpb:
%     spM = penalM(1:(Nt_ts - 1):(allt_lasti), 1:(Nt_ts - 1):(allt_lasti));
    J1fun = zeros(1, Nt + 1);
    % Computation of the weights of the Chebyshev points in integration:
    weightcheb = chebweights(Nt_ts, dt);
    % The weights of the whole grid. The weight of the points on the
    % boundary between 2 time steps is multiplied by 2:
    gridintw = [weightcheb(1:(Nt_ts - 1)), kron(ones(1, Nt - 1), [weightcheb(1)*2, weightcheb(2:(Nt_ts - 1))]), weightcheb(Nt_ts)];
    % The penalty matrix columns are weighted with the weights of integration:
    penalM(:, [1, allt_lasti]) = penalM(:, [1, allt_lasti])*weightcheb(1);
    penalM(:, Nt_ts:(Nt_ts - 1):((Nt - 1)*(Nt_ts - 1) + 1)) = penalM(:, Nt_ts:(Nt_ts - 1):((Nt - 1)*(Nt_ts - 1) + 1))*weightcheb(1)*2;
    % Using the fact that the weights are symmetric around the middle of
    % the domain of integration in any time step:
    for ti = 2:ceil(Nt_ts/2)
        penalM(:, [ti:(Nt_ts - 1):(allt_lasti), (Nt_ts - ti + 1):(Nt_ts - 1):(allt_lasti)]) =...
            penalM(:, [ti:(Nt_ts - 1):(allt_lasti), (Nt_ts - ti + 1):(Nt_ts - 1):(allt_lasti)])*weightcheb(ti);
    end
    maxNiter = 100;
    conv = zeros(1, maxNiter + 1);
    allmniter = 0;
tic
    [allpsi, fieldt(1:Nt), mniter] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
    allmniter = allmniter + mniter;
toc
    fieldt(Nt + 1) = allfield(allt_lasti);
    lastfieldt = fieldt;
    niter = 0;
    getconv();
    relE = tol + 1;
    allfield(allt_lasti) = (-imag(chiT.'*miu*allpsi(:, Nt_ts, Nt))...
        - penalM(allt_lasti, 1:(allt_lasti - 1))*allfield(1:(allt_lasti - 1)))/penalM(allt_lasti, allt_lasti);
    while relE>tol && niter<maxNiter
        tic
        get_chiihterm(allpsi);
        [allchi, allfield(allt_lasti:-1:1), mniter] = solveOCMLSihttp(@ihchifttpMLS1, H0, Edomain, miu, chiT, [T 0],...
            Nt, Nt_ts, Ncheb, tolprop, allfield((allt_lasti):-1:1), allpsi(:, Nt_ts:-1:1, Nt:-1:1),...
            penalM((allt_lasti):-1:1, (allt_lasti):-1:1), chiihterm);
        allmniter = allmniter + mniter;
        [allpsi, allfield, mniter] = solveOCMLSihttp(@ihpsittpMLS1, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop,...
            allfield, allchi(:, Nt_ts:-1:1, Nt:-1:1), penalM);
        allmniter = allmniter + mniter;
        fieldt = allfield(1:(Nt_ts - 1):(allt_lasti)).';
        relE = norm(fieldt - lastfieldt)/norm(fieldt);
        lastfieldt = fieldt;
        niter = niter + 1;
        getconv();
        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    conv = conv(1:(niter+1));
    fieldw = dctI(fieldt);
    niter
    mniterc = allmniter/(2*niter + 1)
    J1 = sum(J1fun);
    grad = zeros(Nt + 1, 1);
    for tsi = 1:Nt
        grad(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*miu*allpsi(:, 1, tsi));
    end
    grad(Nt + 1) = -imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt));
    gradJ2 = -penalM*allfield;
    grad = grad + gradJ2(1:(Nt_ts - 1):(allt_lasti));
    maxgrad = max(grad);
    
    %%% Nested functions: %%%
    
    function wvector = getwvector(winput)
        if length(winput) == 1
            % if winput is a function handle:
            wvector = zeros(1, Nt + 1);
           for wi = 1:(Nt + 1)
                wvector(wi) = winput((wi-1)*dw);
           end
        else
            % if winput is a vector:
            wvector = winput;
        end
    end

    function getconv()
        psi(:, 1:Nt) = allpsi(:, 1, :);
        psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
        evmiut = evmiuMLS(psi, miu);
        evmiuw = dctI(evmiut);
        J1fun = 0.5*evmiuw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
%         % The field vector for integration over time:
%         fieldtint = fieldt*dt;
%         fieldtint([1, Nt + 1]) = fieldtint([1, Nt + 1])/2;
%         conv(:, niter+1) = sum(J1fun) - fieldtint.'*spM*fieldtint;
        conv(:, niter+1) = sum(J1fun) - (allfield.'.*gridintw)*penalM*allfield;
    end

    function get_chiihterm(allpsi)
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)));
        for ti = 1:Nt_ts
            chiihterm(:, ti, 1) = -evmiufil(Nt*(Nt_ts - 1) - ti + 2)*miu*allpsi(:, Nt_ts - ti + 1, Nt);
        end        
        for tsi = 2:Nt
            chiihterm(:, 1, tsi) = chiihterm(:, Nt_ts, tsi - 1);
            for ti = 2:Nt_ts
                chiihterm(:, ti, tsi) = -evmiufil((Nt - tsi + 1)*(Nt_ts - 1) - ti + 2)*miu*allpsi(:, Nt_ts - ti + 1, Nt - tsi + 1);
            end
        end    
    end

end