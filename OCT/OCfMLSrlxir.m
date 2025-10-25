function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1] = OCfMLSrlxir(psi0, E0, Edomain, miu, fguess, ...
    filterE, filtermiu, T, dt, Nt_ts, Ncheb, tol)
%%% The program doesn't work well. I don't know if the reason is the
%%% algorithm, or a bug.
% fguess, filterE, filtermiu, are function handles of w.
    Nt = T/dt;
    dw = pi/T;
%    Nwfield = ceil(maxwE/dw);
%    Nwmiu = ceil(maxwmiu/dw);
    Npsi = length(psi0);
    maxNiter = 20;
    tolprop = 1e-2*tol;
    conv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
    summniterc = 0;
%    allfield = zeros(1, Nt*(Nt_ts - 1) + 1);
    fieldw = zeros(1, Nt + 1);
    evmiut = zeros(1, Nt + 1);
    evmiuw = zeros(1, Nt + 1);
    evmiufil = zeros(1, Nt + 1);
    J1fun = zeros(1, Nt + 1);
    chimiupsi = zeros(1, Nt + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    chiT = zeros(Npsi, 1);
    psi = zeros(Npsi, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
%    vfiltermiu = zeros(1, Nwmiu + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Npsi, Nt_ts, Nt);
    chiihterm = zeros(Npsi, Nt_ts, Nt);
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = 1;
    maxNiter = 20;
    niter = 0;
    allfield = dctIintgrid(fieldw, T, t_ts(1:(Nt_ts-1)));
    [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
    summniterc = summniterc + mniterc;
    getconv(fieldw);
    get_chiihterm(allpsi);
    [allchi, field1, mniterc] = solveOCMLSih(@ihfieldchifMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
        allfield((Nt*(Nt_ts - 1) + 1):-1:1), chiihterm);
    summniterc = summniterc + mniterc;
    nprop = 2;
    chimiupsi = get_cmp(allchi(:, Nt_ts:-1:1, Nt:-1:1), allpsi, miu);
    niter = 1;
    searchfield(@proppsi);
    while relE>tol && niter<maxNiter
%        tic
        niter = niter + 1;
        searchfield(@propchipsi);
        relE = norm(fieldw - lastfieldw)/norm(fieldw);
        lastfieldw = fieldw;
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    niter = niter + 1;
    allfield = dctIintgrid(fieldw, T, t_ts(1:(Nt_ts-1)));
    [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
    summniterc = summniterc + mniterc;
    getconv(fieldw);
    conv = conv(1:(niter+1));
    niter
    weight
    nprop
    fieldt = dctI(fieldw);
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    
    %%% Nested functions: %%%
    function getconv(fieldwint)
        psi(:, 1:Nt) = allpsi(:, 1, :);
        psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
        evmiut = evmiuMLS(psi, miu);
        evmiuw = dctI(evmiut);
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evmiuw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        %integrand = (-0.5*evmiuw.^2.*vfiltermiu,  + J2fun)*dw/(Nt_ts-1);
        %integrand([1, Nt*(Nt_ts - 1)]) = integrand([1, Nt*(Nt_ts-1)])/2;
        %minusJ = sum(integrand); 
        conv(niter + 1) = sum(J1fun) + sum(J2fun);
    end

    function searchfield(propfun)
        flag = 0;
        while flag == 0
            newchimiupsi = propfun();
            newfieldw = dctI(newchimiupsi).*vfilterE;
            tryfieldw = (1 - weight)*fieldw + weight*newfieldw;
            getconv(tryfieldw);
            if conv(niter + 1) > conv(niter)
                flag = 1;                
            else
                weight = weight/2;
            end
        end    
        fieldw = tryfieldw;
        chimiupsi = newchimiupsi;
    end

    function newchimiupsi = proppsi()
        [allpsi, newchimiupsi, mniterc] = solveOCfMLS(@ihpsifMLSrlx, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, ...
            chimiupsi, allchi(:, Nt_ts:-1:1, Nt:-1:1), vfilterE, fieldw, weight);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
    end

    function newchimiupsi = propchipsi()
        get_chiihterm(allpsi);
        [allchi, newchimiupsi, mniterc] = solveOCfMLS(@ihchifMLSrlx, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop, ...
            chimiupsi, allpsi(:, Nt_ts:-1:1, Nt:-1:1), vfilterE, chiihterm, fieldw, weight);
        summniterc = summniterc + mniterc;
        [allpsi, newchimiupsi, mniterc] = solveOCfMLS(@ihpsifMLSrlx, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, ...
            newchimiupsi, allchi(:, Nt_ts:-1:1, Nt:-1:1), vfilterE, fieldw, weight);
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
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

function chimiupsi = get_cmp(chi, psi, miu)
    [Npsi, Nt_ts, Nt] = size(psi);
    chimiupsi = zeros(1, Nt + 1);
     for tsi = 1:Nt
         chimiupsi(tsi) = -imag(chi(:, 1, tsi)'*miu*psi(:, 1, tsi));
     end
     chimiupsi(Nt+1) = -imag(chi(:, Nt_ts, Nt)'*miu*psi(:, Nt_ts, Nt));
end