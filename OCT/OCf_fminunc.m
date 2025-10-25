function [fieldt, fieldw, psi, evmiut, evmiuw, conv, niter, mallniterc, J1, output] = OCf_fminunc(psi0, H0, Edomain, miu,...
    fguess, filterE, filtermiu, T, dt, Nt_ts, Ncheb, tol, maxNiter)
% It works for rectangular filtration of the field only.
    Nt = T/dt;
    dw = pi/T;
    dctfactor = T/(sqrt(Nt*pi));
    dim = length(psi0);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    tolprop = 1e-3*tol;
    summniterc = 0;
    allfield = zeros(allt_lasti, 1);
    allevmiut = zeros(allt_lasti, 1);
    evmiufil = zeros(allt_lasti + 1, 1);
    evmiuw = zeros(1, Nt + 1);
    J1fun = zeros(Nt + 1, 1);
    chimiupsi = zeros(allt_lasti, 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
%     allt_grid = [kron((0:dt:(T - dt)).', ones(Nt_ts - 1, 1)) + kron(ones(Nt, 1), t_ts(1:(Nt_ts - 1))); T];
    w = (0:pi/T:pi/dt).';
    vfilterE = filterE(w);
%    [iEnz, ~, vfilterEnz] = find(vfilterE);
    iEnz = find(vfilterE>=tolprop*max(vfilterE)/10);
    vfilterEnz = vfilterE(iEnz);
    Nwnz = length(iEnz);
    fieldw = zeros(Nt + 1, 1);
    if length(fguess) == 1
        % if fguess is a function handle:
        for wnzi = iEnz
            fieldw(wnzi) = fguess(w(wnzi));
        end
    else
        % if fguess is a vector:
        fieldw(iEnz) = fguess(iEnz);
    end
    fieldwnz = fieldw(iEnz);
    vfiltermiu = filtermiu(w);
    allpsi = zeros(dim, allt_lasti);
    allmiupsi = zeros(dim, allt_lasti);
    chiT = zeros(dim, 1);
    chiihterm = zeros(dim, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    igweights = integw(Nt_ts:2*(Nt_ts - 1))/dt;
%    integwnz = integw(iEnz);
    nprop = 0;
    options = optimoptions('fminunc', 'algorithm', 'quasi-newton', 'OptimalityTolerance', tol, 'StepTolerance', tol, 'SpecifyObjectiveGradient',true,...
        'OutputFcn', @getconv, 'MaxIter', maxNiter);
    figure
    [fieldwnz, ~, ~, output] = fminunc(@Jeval, fieldwnz, options);
    output
    nprop
    fieldw = fieldw.';
%    maxgrad = max(abs(minusgrad));
    %conv = -minus_conv;
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = allfield(1:(Nt_ts - 1):allt_lasti).';
    evmiut = allevmiut(1:(Nt_ts - 1):allt_lasti).';
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    beep
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad] = Jeval(fieldwnz)
        fieldw(iEnz) = fieldwnz;
        allfield = dctIintgrid(fieldw, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield.');
        summniterc = summniterc + mniterc;
        allmiupsi = miu*allpsi;
        for allti = 1:allt_lasti
% The values of the expectation value are supposed to be real:
            allevmiut(allti) = real(allpsi(:, allti)'*allmiupsi(:, allti));
        end
        evmiuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        get_chiihterm(allmiupsi);
        [allchi, ~, mniterc] = solveOCMLSih1(@ihalltchiMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1), chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        J1fun = 0.5*evmiuw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwnz.^2./vfilterEnz*dw;
        if iEnz(1) == 1
            J2fun(1) = J2fun(1)/2;
        end
        if iEnz(Nwnz) == (Nt + 1)
            J2fun(Nwnz) = J2fun(Nwnz)/2;
        end
        %J2fun(iEnz==1 || iEnz==(Nt + 1)) = J2fun(iEnz==1 || iEnz==(Nt + 1))/2;
        minusJ = -(sum(J1fun) + sum(J2fun));
        for allti = 1:allt_lasti
            chimiupsi(allti) = imag(allchi(:, allt_lasti - allti + 1)'*(miu*allpsi(:, allti)));
        end
        dct_chimiupsi = dctfactor*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights);
        minusgrad = 2*dw*(fieldwnz./vfilterEnz + dct_chimiupsi(iEnz));
        if iEnz(1) == 1
            minusgrad(1) = minusgrad(1)/2;
        end
        if iEnz(Nwnz) == (Nt + 1)
            minusgrad(Nwnz) = minusgrad(Nwnz)/2;
        end
    end


    function get_chiihterm(allmiupsi)
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        chiihterm = -allmiupsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1), 0, allt_lasti, allt_lasti);     
    end

    function stop = getconv(fieldwnz, optimValues, state)
        stop = false;
        niter = optimValues.iteration;
        conv(optimValues.iteration + 1) = -optimValues.fval;        
        plot(0:niter, conv)
        drawnow
    end

end