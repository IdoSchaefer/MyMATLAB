function [fieldt, fieldw, psi, evmiut, evmiuw, evmiutcut, evmiuwcut, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfMLSbd0(psi0, E0, Edomain, miu,...
    fguess, filterE, filtermiu, tcutfun, Dtcutfun, iweight, T, dt, Nt_ts, Ncheb, tol)
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
% tol: the tolerance of the optimization process.
    Nt = T/dt;
    dw = pi/T;
    w = 0:dw:pi/dt;
    dctfactor = T/(sqrt(Nt*pi));
%    Nwfield = ceil(maxwE/dw);
%    Nwmiu = ceil(maxwmiu/dw);
    Npsi = length(psi0);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
    summniterc = 0;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    allfield = zeros(1, allt_lasti);
    if length(fguess) == 1
        % if fguess is a function handle:
        fieldw = zeros(1, Nt + 1);
        for wi = 1:(Nt + 1)
            fieldw(wi) = fguess(w(wi));
        end
    else
        % if fguess is a vector:
        fieldw = fguess;
    end
%    evmiut = zeros(1, Nt + 1);
    evmiuallt = zeros(1, allt_lasti);
    evmiualltcut = zeros(1, allt_lasti);
    evmiuwcut = zeros(1, Nt + 1);
    evmiufil = zeros(1, Nt + 1);
    J1fun = zeros(1, Nt + 1);
    chimiupsi = zeros(1, allt_lasti);
    allmiupsi = zeros(1, allt_lasti);
    t = 0:dt:T;
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt = [kron(t(1:Nt), ones(1, Nt_ts - 1)), T] + [kron(ones(1, Nt), t_ts(1:(Nt_ts - 1))), 0];
    tcutallt = tcutfun(allt);
    Dtcutallt = Dtcutfun(allt);
%     integw = chebgridw(Nt, Nt_ts, dt);
%     igweights = integw(Nt_ts:2*(Nt_ts - 1))/dt;    
%    igweights = chebweights(Nt_ts, 1);
    alligweights = chebweights_gen(Nt_ts, [0 dt], t_ts(2:Nt_ts));
    igweights = [alligweights(1, Nt_ts - 1)*2, alligweights(2:(Nt_ts - 1), Nt_ts - 1).']/dt;
    chiT = zeros(Npsi, 1);
    vfilterE = zeros(1, Nt + 1);
%    vfiltermiu = zeros(1, Nwmiu + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Npsi, allt_lasti);
    chiihterm = zeros(Npsi, allt_lasti);
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = iweight;
    niter = 0;
    getpsi(fieldw);
    nprop = 1;
    while relE>tol && niter<maxNiter
%        tic
        get_chiihterm(allmiupsi);
        [allchi, field1, mniterc] = solveOCMLSih1(@ihalltchiMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield((Nt*(Nt_ts - 1) + 1):-1:1), chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        for allti = 1:allt_lasti
            chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*allmiupsi(:, allti));
        end
        niter = niter + 1;
        flag = 0;
        newfieldw = dctfactor*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights).*vfilterE;
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
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = dctI(fieldw)/dctfactor;
    evmiut = evmiuallt(1:(Nt_ts - 1):allt_lasti);
    evmiuw = dctIfrom_ig(evmiuallt, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;    
    evmiutcut = evmiualltcut(1:(Nt_ts - 1):allt_lasti);
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    for allti = 1:allt_lasti
        chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*allmiupsi(:, allti));
    end
    grad = zeros(1, Nt + 1);
    grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
    grad = grad + 2*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    maxgrad = max(abs(grad));
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        allmiupsi = miu*allpsi;
        for allti = 1:allt_lasti
            evmiuallt(allti) = allpsi(:, allti)'*allmiupsi(:, allti);
        end
        integDtcut = 0;
        evmiualltcut(1) = evmiuallt(1)*tcutallt(1);
        for tsi = 1:Nt
            allti_ts = (tsi - 1)*(Nt_ts - 1) + (1:Nt_ts);
            integrand_ts = evmiuallt(allti_ts).*Dtcutallt(allti_ts);
            ts_integral = integrand_ts*alligweights;
            evmiualltcut(allti_ts(2:Nt_ts)) = evmiuallt(allti_ts(2:Nt_ts)).*tcutallt(allti_ts(2:Nt_ts)) - integDtcut - ts_integral;
            integDtcut = integDtcut + ts_integral(Nt_ts - 1);
        end
        evmiuwcut = dctIfrom_ig(evmiualltcut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evmiuwcut.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        %integrand = (-0.5*evmiuw.^2.*vfiltermiu,  + J2fun)*dw/(Nt_ts-1);
        %integrand([1, Nt*(Nt_ts - 1)]) = integrand([1, Nt*(Nt_ts-1)])/2;
        %minusJ = sum(integrand); 
        conv(niter + 1) = sum(J1fun) + sum(J2fun);
    end

    function get_chiihterm(allmiupsi)
        evmiufilw = evmiuwcut.*vfiltermiu;
        evmiufil = dctIintgrid(evmiufilw, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        sinc_transform = sinct_intgrid(evmiufilw, T, t_ts(1:(Nt_ts-1)));
        chiihterm = -(allmiupsi(:, allt_lasti:-1:1)*spdiags((evmiufil(allt_lasti:-1:1).*tcutallt(allt_lasti:-1:1) ...
            + (sinc_transform(allt_lasti:-1:1) - sinc_transform(allt_lasti)).*Dtcutallt(allt_lasti:-1:1)).', 0, allt_lasti, allt_lasti)); 
    end

end