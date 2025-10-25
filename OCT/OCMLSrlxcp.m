function [field, psi, relE, conv, niter, convprop, npropconv, mallniterc, J1, maxgrad, weight] = OCMLSrlxcp(psi0, target, E0, Edomain,...
    miu, fguess, Epenal, T, dt, Nt_ts, Ncheb, tol)
% The program uses a relxation process for an optimal control of a state to
% state transition.
% Returns also convprop - the convergens vs. the number of propagations -
% in npropconv.
% fguess is a function handle of t.
    Nt = T/dt;
    Npsi = length(psi0);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    convprop = zeros(1, maxNiter + 1);
    npropconv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
    summniterc = 0;
    J1 = 0;
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    allfield = zeros(1, allt_lasti);
    newallfield = zeros(1, allt_lasti);
    for tsi = 1:Nt
        for ti = 1:(Nt_ts - 1)
            allfield((tsi - 1)*(Nt_ts - 1) + ti) = fguess((tsi - 1)*dt + t_ts(ti));
        end
    end
    allfield(allt_lasti) = fguess(T);
    lastfield = allfield(1:(Nt_ts - 1):allt_lasti);
    psi = zeros(Npsi, Nt + 1);
    allpsi = zeros(Npsi, Nt_ts, Nt);
    integw = chebgridw(Nt, Nt_ts, dt);
    overlap = 0;
    relE = tol + 1;
    weight = 1;
    niter = 0;
    getpsi(allfield);
    nprop = 1;
    convprop(niter + 1) = conv(niter + 1);
    npropconv(niter + 1) = 1;
    while relE>tol && niter<maxNiter
%        tic
        chiT = overlap*target;
        [allchi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1));
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        niter = niter + 1;
        flag = 0;
        for tsi = 1:Nt
            for ti = 1:(Nt_ts - 1)
                newallfield((tsi - 1)*(Nt_ts - 1) + ti) = -imag(allchi(:, Nt_ts - ti + 1, Nt - tsi + 1)'*miu*allpsi(:, ti, tsi))/Epenal;
            end
        end
        newallfield(allt_lasti) = -imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt))/Epenal;
        while flag == 0
            tryallfield = (1 - weight)*allfield + weight*newallfield;
            getpsi(tryallfield);
            nprop = nprop + 1;
            if conv(niter + 1) > conv(niter)
                flag = 1;                
            else
                weight = weight/2;
            end
        end
        allfield = tryallfield;
        field = allfield(1:(Nt_ts - 1):allt_lasti);
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        convprop(niter + 1) = conv(niter + 1);
        npropconv(niter + 1) = nprop;
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    conv = conv(1:(niter+1));
    convprop = convprop(1:(niter+1));
    npropconv = npropconv(1:(niter+1));
    figure
    plot(npropconv, convprop)
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    niter
    weight
    nprop
    mallniterc = summniterc/nprop;
    grad = 2*(newallfield - allfield)*Epenal; 
    maxgrad = max(abs(grad));
    
    %%% Nested functions: %%%
    function getpsi(allfieldint)
        [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfieldint);
        summniterc = summniterc + mniterc;
        psiT = allpsi(:, Nt_ts, Nt);
        overlap = target'*psiT;
        J1 = overlap*conj(overlap);
        conv(niter + 1) = J1 - Epenal*sum(integw.*allfieldint.^2);
    end

end