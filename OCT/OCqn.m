function [allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn(psi0, target, H0, Edomain, miu,...
    fguess, Epenal, options, T, dt, Nt_ts, Ncheb, tol, maxNiter)
    Nt = round(T/dt);
    dim = length(psi0);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    tolprop = 1e-3*tol;
    summniterc = 0;
    allfield = zeros(allt_lasti, 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    if length(fguess) == 1
        % if fguess is a function handle:
        for tsi = 1:Nt
            for ti = 1:(Nt_ts - 1)
                allfield((tsi - 1)*(Nt_ts - 1) + ti) = fguess((tsi - 1)*dt + t_ts(ti));
            end
        end
        allfield(Nt*(Nt_ts - 1) + 1) = fguess(Nt*dt);
    else
        % if fguess is a vector:
        allfield = fguess;
    end
    J1 = 0;
    allpsi = zeros(dim, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    nprop = 0;
    if isempty(options)
        options = optionsOCqn(tol, maxNiter);
    end
    if isempty(options.invHess0)
        options.invHess0 = diag(1./(2*Epenal*integw));
    end
    %options.sigma = 0.1;
    %options.tau1 = 2;
    [allfield, ~, minus_grad, niter, ~, ~, dif_field, minus_conv, alpha, invHess] = quasiNewton(@Jeval, allfield, options);
    allfield = allfield.';
    field = allfield(1:(Nt_ts - 1):allt_lasti);
    %grad = -minus_grad;
    maxgrad = max(abs(minus_grad));
    conv = -minus_conv;
    relE = norm(dif_field)/norm(allfield);
    niter
    nprop
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    mallniterc = summniterc/nprop;
    beep
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad] = Jeval(allfield)
        [allpsi, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        overlap = target'*allpsi(:, allt_lasti);
        chiT = overlap*target;
        [allchi, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1));
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        J1 = overlap*conj(overlap);
        minusJ = -J1 + Epenal*sum(integw.'.*allfield.^2);
        minusgrad = zeros(allt_lasti, 1);
        for allti = 1:allt_lasti
            minusgrad(allti) = 2*integw(allti)*(allfield(allti)*Epenal + imag(allchi(:, allt_lasti - allti + 1)'*(miu*allpsi(:, allti))));
        end
    end

end