function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad] = OCfMLSaloud(psi0, E0, Edomain, miu, ...
    aloudwf, forbwf, fguess, filterE, filtermiu, penalforb, iweight, T, dt, Nt_ts, Ncheb, tol)
% The program finds an optimal field for frequency control. It's intended
% for a many level system.
% It uses a relaxation process.
% psi0:  The initial state
% E0: The eigenenergies of H0
% miu: The dipole operator matrix
% fguess, filterE, filtermiu, are function handles of w (the angular freqency), 
% or row vectors, when the w values are: 0:pi/T:pi/dt
% fguess: the initial guess for the field.
% filterE: the function that filters the field in the frequncy domain,
% multiplied by 1/(penalty factor)
% filtermiu: the function that filters the dipole eigenvalues in the
% frequncy domain.
% penalforb: a function handle of the number of column in forbwf. Defines
% the penalty on the forbiden wave functions - each function may have a different penalty.
% tol: the tolerance of the optimization process.
    Nt = T/dt;
    dw = pi/T;
    dctfactor = T/(sqrt(Nt*pi));
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    Npsi = length(psi0);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
    projaloud = aloudwf*aloudwf';
    miualoud = projaloud*miu*projaloud;
%     projforb = eye(Npsi) - projaloud;
    Nforb = size(forbwf, 2);
    penalforbv = penalforb(1:Nforb);
    projforb = forbwf*diag(penalforbv)*forbwf';
    allmiualpsi = zeros(Npsi, allt_lasti);
    miualpsi = zeros(Npsi, Nt + 1);
    allpfpsi = zeros(Npsi, allt_lasti);
    allevprojforb = zeros(1, allt_lasti);
    summniterc = 0;
    allfield = zeros(1, Nt*(Nt_ts - 1) + 1);
    if length(fguess) == 1
        % if fguess is a function handle:
        fieldw = zeros(1, Nt + 1);
        for wi = 1:(Nt + 1)
            fieldw(wi) = fguess((wi-1)*dw);
        end
    else
        % if fguess is a vector:
        fieldw = fguess;
    end
    evmiualt = zeros(1, Nt + 1);
    evmiualw = zeros(1, Nt + 1);
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
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Npsi, Nt_ts, Nt);
    chiihterm = zeros(Npsi, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = iweight;
    niter = 0;
    getpsi(fieldw);
    nprop = 1;
    while relE>tol && niter<maxNiter
%        tic
        get_chiihterm(allpsi);
        %chiT = allpsi(:, Nt_ts, Nt);
        [allchi, field1, mniterc] = solveOCMLSih(@ihalltchiMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield((Nt*(Nt_ts - 1) + 1):-1:1), chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        for tsi = 1:Nt
            chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*miu*allpsi(:, 1, tsi));
        end
        chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt));
        niter = niter + 1;
        flag = 0;
        newfieldw = dctfactor*dctI(chimiupsi).*vfilterE;
        while flag == 0
            tryfieldw = (1 - weight)*fieldw + weight*newfieldw;
            getpsi(tryfieldw);
            nprop = nprop + 1;
            if conv(niter + 1) > conv(niter)
                flag = 1;                
            else
                weight = weight/2;
            end
        end
        fieldw = tryfieldw;
        relE = norm(fieldw - lastfieldw)/norm(fieldw);
        lastfieldw = fieldw;
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    conv = conv(1:(niter+1));
    niter
    weight
    nprop
    fieldt = dctI(fieldw)/dctfactor;
    evmiut = evmiuMLS(psi, miu);
    evmiuw = dctI(evmiut)*dctfactor;
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    for tsi = 1:Nt
        chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*miu*allpsi(:, 1, tsi));
    end
    chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt));
    grad = zeros(1, Nt + 1);
    grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
    grad = grad + 2*dctI(chimiupsi)*dctfactor;
    maxgrad = max(abs(grad));
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        psi(:, 1:Nt) = allpsi(:, 1, :);
        psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
%        evmiualt = evmiuMLS(psi, miualoud);
        miualpsi = miualoud*psi;
        for tsi = 1:(Nt + 1)
            evmiualt(tsi) = psi(:, tsi)'*miualpsi(:, tsi);
        end
        evmiualw = dctI(evmiualt)*dctfactor;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evmiualw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
%         Jforbfun = -evmiuMLS(psi, projforb)*dt;
%         Jforbfun([1, Nt + 1]) = Jforbfun([1, Nt + 1])/2;
        for ti = 1:(Nt_ts - 1)
            psiti(:, :) = allpsi(:, ti, :);
            allpfpsi(:, ti:(Nt_ts - 1):(allt_lasti - 1)) = projforb*psiti;
            for tsi = 1:Nt
                allevprojforb((tsi - 1)*(Nt_ts - 1) + ti) = allpsi(:, ti, tsi)'*allpfpsi(:, (tsi - 1)*(Nt_ts - 1) + ti);
            end
        end
        allpfpsi(:, allt_lasti) = projforb*allpsi(:, Nt_ts, Nt);
        allevprojforb(allt_lasti) = allpsi(:, Nt_ts, Nt)'*allpfpsi(:, allt_lasti);
        conv(niter + 1) = sum(J1fun) + sum(J2fun) - sum(integw.*allevprojforb);
    end

    function get_chiihterm(allpsi)
        evmiufil = dctIintgrid(evmiualw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        allmiualpsi(:, 1:(Nt_ts - 1):allt_lasti) = miualpsi;
        for ti = 2:(Nt_ts - 1)
            psiti(:, :) = allpsi(:, ti, :);
            allmiualpsi(:, ti:(Nt_ts - 1):(allt_lasti - 1)) = miualoud*psiti;
        end
        allmiualpsi(:, allt_lasti) = miualoud*allpsi(:, Nt_ts, Nt);
        chiihterm = allpfpsi(:, allt_lasti:-1:1) -...
            allmiualpsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1).', 0, allt_lasti, allt_lasti); 
%         for ti = 1:Nt_ts
%             chiihterm(:, ti, 1) = (projforb - evmiufil(Nt*(Nt_ts - 1) - ti + 2)*miualoud)*allpsi(:, Nt_ts - ti + 1, Nt);
%         end
%         for tsi = 2:Nt
%             chiihterm(:, 1, tsi) = chiihterm(:, Nt_ts, tsi - 1);
%             for ti = 2:Nt_ts
%                 chiihterm(:, ti, tsi) = (projforb - evmiufil((Nt - tsi + 1)*(Nt_ts - 1) - ti + 2)*miualoud)...
%                     *allpsi(:, Nt_ts - ti + 1, Nt - tsi + 1);
%             end
%         end    
    end

end