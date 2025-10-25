function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(psi0, Vf, m, Edomain,...
    xdomain, miufunx, aloudwf, forbwf, fguess, filterE, filtermiu, orthpenal, penalforb, iweight, T, dt, Nt_ts, Ncheb, tol)
% The program finds an optimal field for frequency control. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% It uses a relaxation process.
%%%% There are mistakes in the computation of |chi(T)>:
%%%% [H(T), projaloud*miux*projaloud] has to be used instead of [K, miux]. Not very important.%%%%%%%
% psi0:  The initial state
% Vf: the potential energy of H0, a function handle of the form: @(x)
% m: The mass
% xdomain: The domain of the equally spaced grid of x: [min_x max_x)
% miufunx: a function of x that represents miu.
% fguess, filterE, filtermiu, are function handles of w (the angular freqency), 
% or row vectors, when the w values are: 0:pi/T:pi/dt
% fguess: the initial guess for the field.
% filterE: the function that filters the field in the frequncy domain,
% multiplied by 1/(penalty factor)
% filtermiu: the function that filters the dipole eigenvalues in the
% frequncy domain.
% tol: the tolerance of the optimization process.
    Nt = T/dt;
    dw = pi/T;
    dctfactor = T/(sqrt(Nt*pi));
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    Nx = length(psi0);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/(2*m);
    miux = miufunx(x);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    projaloud = aloudwf*aloudwf';
    miualoud = projaloud*spdiags(miux, 0, Nx, Nx)*projaloud;
%    projforb = eye(Nx) - projaloud;
    Nforb = size(forbwf, 2);
    penalforbv = penalforb(1:Nforb);
    projforb = forbwf*diag(penalforbv)*forbwf';
    allmiualpsi = zeros(Nx, allt_lasti);
    miualpsi = zeros(Nx, Nt + 1);
    allpfpsi = zeros(Nx, allt_lasti);
    allevprojforb = zeros(1, allt_lasti);
%    evpT = 0;
    evcomT = 0;
    compsiT = zeros(Nx, 1);
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
    psi = zeros(Nx, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
%    vfiltermiu = zeros(1, Nwmiu + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Nx, Nt_ts, Nt);
    chiihterm = zeros(Nx, allt_lasti);
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
%        chiT = -ppenal*evpT*ifft(p.*fft(allpsi(:, Nt_ts, Nt)))/m^2;
        chiT = orthpenal*evcomT*compsiT;
        [allchi, field1, mniterc] = solveOCih(@ihalltchimiux, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield((Nt*(Nt_ts - 1) + 1):-1:1), miux, chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        for tsi = 1:Nt
            chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*(miux.*allpsi(:, 1, tsi)));
        end
        chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*(miux.*allpsi(:, Nt_ts, Nt)));
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
    evmiut = evmiu(psi, miux);
    evmiuw = dctI(evmiut)*dctfactor;
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    for tsi = 1:Nt
        chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*(miux.*allpsi(:, 1, tsi)));
    end
    chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*(miux.*allpsi(:, Nt_ts, Nt)));
    grad = zeros(1, Nt + 1);
    grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
    grad = grad + 2*dctI(chimiupsi)*dctfactor;
    maxgrad = max(abs(grad));
%    evp = evpT
    evorth = 1i*evcomT
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCih(@ihfieldmiux, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield, miux);
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
%        allevprojforb = -evallpsi(allpsi, projforb);
        for ti = 1:(Nt_ts - 1)
            psiti(:, :) = allpsi(:, ti, :);
            allpfpsi(:, ti:(Nt_ts - 1):(allt_lasti - 1)) = projforb*psiti;
            for tsi = 1:Nt
                allevprojforb((tsi - 1)*(Nt_ts - 1) + ti) = allpsi(:, ti, tsi)'*allpfpsi(:, (tsi - 1)*(Nt_ts - 1) + ti);
            end
        end
        allpfpsi(:, allt_lasti) = projforb*allpsi(:, Nt_ts, Nt);
        allevprojforb(allt_lasti) = allpsi(:, Nt_ts, Nt)'*allpfpsi(:, allt_lasti);
        psiT = psi(:, Nt + 1);
        compsiT = ifft(K.*fft(miux.*psiT)) - miux.*ifft(K.*fft(psiT));
        evcomT = psiT'*compsiT;
%        conv(niter + 1) = sum(J1fun) + sum(J2fun) - sum(integw.*allevprojforb) - 0.5*ppenal*(evpT/m)^2;
        conv(niter + 1) = sum(J1fun) + sum(J2fun) - sum(integw.*allevprojforb) + 0.5*orthpenal*evcomT^2;
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