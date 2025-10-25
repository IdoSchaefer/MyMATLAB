function [allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0, target, H0, Edomain, miu,...
    fguess, filterE, options, T, dt, Nt_ts, Ncheb, tol, maxNiter)
    Nt = T/dt;
    dim = length(psi0);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    tolprop = 1e-3*tol;
    summniterc = 0;
    tcheb = -cos(((1:Nt_ts).' - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt_grid = [kron((0:dt:(T - dt)).', ones(Nt_ts - 1, 1)) + kron(ones(Nt, 1), t_ts(1:(Nt_ts - 1))); T]; 
    vfilterE = filterE(allt_grid);
%    [iEnz, ~, vfilterEnz] = find(vfilterE);
    iEnz = find(vfilterE>=tolprop*max(vfilterE)/10);
    vfilterEnz = vfilterE(iEnz);
    Nalltnz = length(iEnz);
    %allt_grid_nz = allt_grid(iEnz);
    allfield = zeros(allt_lasti, 1);
    if length(fguess) == 1
        % if fguess is a function handle:
        for alltnzi = iEnz
            allfield(alltnzi) = fguess(allt_grid(alltnzi));
        end
    else
        % if fguess is a vector:
        allfield(iEnz) = fguess(iEnz);
    end
    allfieldnz = allfield(iEnz);
    %allfield = zeros(allt_lasti, 1);
    %allfield(iEnz) = allfieldnz;
    J1 = 0;
    allpsi = zeros(dim, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    integwnz = integw(iEnz).';
    nprop = 0;
    if isempty(options)
        options = optionsOCqn(tol, maxNiter);
    end
    if isempty(options.invHess0)
        options.invHess0 = diag(vfilterEnz./(2*integwnz));
    end
    %options.sigma = 0.1;
    %options.tau1 = 2;
    [~, ~, minus_grad, niter, ~, ~, dif_field, minus_conv, alpha, invHess] = quasiNewton(@Jeval, allfieldnz, options);
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
    
    function [minusJ, minusgrad] = Jeval(allfieldnz)
        allfield(iEnz) = allfieldnz;
        [allpsi, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        overlap = target'*allpsi(:, allt_lasti);
        chiT = overlap*target;
        [allchi, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1));
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        J1 = overlap*conj(overlap);
        minusJ = -J1 + sum(integwnz.*allfieldnz.^2./vfilterEnz);
        minusgrad = zeros(Nalltnz, 1);
        for gradi = 1:Nalltnz
            minusgrad(gradi) = 2*integwnz(gradi)*(allfieldnz(gradi)./vfilterEnz(gradi) +...
                imag(allchi(:, allt_lasti - iEnz(gradi) + 1)'*(miu*allpsi(:, iEnz(gradi)))));
        end
    end

end