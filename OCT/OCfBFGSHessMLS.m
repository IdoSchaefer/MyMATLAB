function [fieldt, fieldw, psi, evmiut, evmiuw, conv, niter, mallniterc] = OCfBFGSHessMLS(psi0, E0, Edomain, miu, fguess, filterE, filtermiu, ...
    maxwE, T, dt, Nt_ts, Ncheb, tol)
% fguess, filterE, filtermiu, are function handles of w.
% fguess can be a vector.
    Nt = T/dt;
    dw = pi/T;
    dctfactor = T/(sqrt(Nt*pi));
    Nwfield = ceil(maxwE/dw);
%    Nwmiu = ceil(maxwmiu/dw);
    Npsi = length(psi0);
    maxNiter = 400;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
    summniterc = 0;
    nprop = 0;
    if length(fguess) == 1
        % if fguess is a function handle:
        fieldw = zeros(1, Nwfield + 1);
        for wi = 1:(Nwfield + 1)
            fieldw(wi) = fguess((wi-1)*dw);
        end
    else
        % if fguess is a vector:
        fieldw = fguess;
    end
    evmiut = zeros(1, Nt + 1);
    evmiuw = zeros(1, Nt + 1);
    chimiupsi = zeros(1, Nt + 1);
%    allfieldt = zeros(1, Nt*(Nt_ts - 1) + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    chiT = zeros(Npsi, 1);
    psi = zeros(Npsi, Nt + 1);
    vfilterE = zeros(1, Nwfield + 1);
%    vfiltermiu = zeros(1, Nwmiu + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nwfield + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    for wi = 1:(Nt + 1) %(Nwmiu + 1)
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
    end
    chiihterm = zeros(Npsi, Nt_ts, Nt);
%    fieldt = dctI(fieldw);
    options = optimset('gradObj', 'on', 'Hessian', 'on', 'OutputFcn', @getconv, 'TolX', tol, 'MaxIter', maxNiter);  
%    options = optimset('gradObj', 'on', 'OutputFcn', @getconv, 'MaxIter', maxNiter);  
    [fieldw, lastJval, exitflag, output] = fminunc(@Jeval, fieldw, options);
    %exitflag
    output
    nprop
    niter = output.iterations;
    conv = conv(1:niter+1);
%     allfield = dctIgrid([fieldw, zeros(1, Nt - Nwfield)], T, t_ts(1:(Nt_ts-1)));
%     allpsi = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
%     psi(:, 1:Nt) = allpsi(:, 1, :);
%     psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
%     evmiut = evmiuMLS(psi, miu);
%    fieldt(Nt + 1) = allfield(Nt*(Nt_ts - 1) + 1);
    fieldt = dctI([fieldw, zeros(1, Nt - Nwfield)])/dctfactor;
    mallniterc = summniterc/nprop;
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad, minusHess] = Jeval(fieldw)
        allfield = dctIintgrid([fieldw, zeros(1, Nt - Nwfield)], T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        get_chiihterm(allpsi);
        [allchi, field1, mniterc] = solveOCMLSih(@ihfieldchifMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield((Nt*(Nt_ts - 1) + 1):-1:1), chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        notsmall = (abs(fieldw)>10*eps) + (abs(vfilterE)>10*eps);
%        J2fun = interpft(fieldw./vfilterE, Nt*(Nt_ts-1));
        J1fun = real(-0.5*evmiuw.^2.*vfiltermiu*dw);
%        J1fun([1, Nwmiu + 1]) = J1fun([1, Nwmiu + 1])/2;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
%         J2fun = fieldw(notsmall == 2).^2./vfilterE(notsmall == 2)*dw;
%         J2fun([1, end]) = J2fun([1, end])/2;
        J2fun = real(fieldw.^2./vfilterE*dw);
        J2fun([1, Nwfield + 1]) = J2fun([1, Nwfield + 1])/2;
        J2fun = J2fun(notsmall == 2);
        %integrand = (-0.5*evmiuw.^2.*vfiltermiu,  + J2fun)*dw/(Nt_ts-1);
        %integrand([1, Nt*(Nt_ts - 1)]) = integrand([1, Nt*(Nt_ts-1)])/2;
        %minusJ = sum(integrand); 
        minusJ = sum(J1fun) + sum(J2fun);
        minusgrad = zeros(1, Nwfield + 1);
        for tsi = 1:Nt
            chimiupsi(tsi) = imag(allchi(:, Nt_ts, Nt - tsi + 1)'*miu*allpsi(:, 1, tsi));
        end
        chimiupsi(Nt + 1) = imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt));
%        minusgrad(notsmall == 2) = 2*fieldw./vfilterE;
%        minusgrad(notsmall ~= 2) = 0
        minusgrad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
        chimiupsiw = dctI(chimiupsi)*dctfactor;
        minusgrad = minusgrad + 2*chimiupsiw(1:(Nwfield + 1));
        % An approximation of the Hessian:
        minusHess = diag(2./vfilterE);
    end
    
    function get_chiihterm(allpsi)
        psi(:, 1:Nt) = allpsi(:, 1, :);
        psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
        evmiut = evmiuMLS(psi, miu);
        evmiuw = dctI(evmiut)*dctfactor;
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
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
    
    function stop = getconv(fieldw, optimValues, state)
        stop = false;
        optimValues.iteration
        conv(optimValues.iteration + 1) = -optimValues.fval;
    end
end