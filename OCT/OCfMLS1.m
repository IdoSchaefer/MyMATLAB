function [fieldw, psi, evmiut, evmiuw, relE, conv, niter, J1] = OCfMLS1(psi0, H0, miumat, fguess, filterE, filtermiu, T, dt, tol)
% fguess: a function handle. A function of the frequency.
    Nt = T/dt;
    dw = pi/T;
    Npsi = length(psi0);
%    W = pi/dt;
    fieldw = zeros(1, Nt + 1);
%    fieldt = zeros(1, Nt + 1);
    chimiupsi = zeros(1, Nt + 1);
%    lastfieldw = zeros(1, Nt + 1);
    chi = zeros(Npsi, Nt + 1);
%    chi0 = zeros(Npsi, 1);
    chiT = zeros(Npsi, 1);
%    chi0 = psi0;
    psi = zeros(Npsi, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
        fieldw(wi) = fguess((wi-1)*dw);
    end
    fieldt = dctI(fieldw);
    lastfieldw = fieldw;
    evmiut = zeros(1, Nt + 1);
    evmiuw = zeros(1, Nt + 1);
    evmiufil = zeros(1, Nt + 1);
    relE = tol + 1;
    conv = zeros(1, 1001);
    weight = 1;
    maxNiter = 20;
    niter = 0;
    getpsi(fieldw);
    prop = 1;
    while relE>tol && niter<maxNiter
%        tic
%        fieldw = dct(chimiupsi).*vfilterE;
        [tchi, chisol] = ode45(@OCfMLSchi, [T 0], chiT, [], H0, miumat, fieldt, psi, evmiufil, Npsi, Nt, dt);
        prop = prop + 1;
        chisol = chisol.';
        chisol = chisol(:, end:-1:1);
        tchi = tchi(end:-1:1);
        for ti = 1:(Nt+1)
            chi(:, ti) = RK4interp(chisol, tchi, (ti - 1)*dt);
        end
        for ti = 1:Nt
            chimiupsi(ti) = -imag(chi(:, ti)'*miumat*psi(:, ti));
        end
        niter = niter + 1;
        flag = 0;
        newfieldw = dctI(chimiupsi).*vfilterE;
        while flag == 0
            tryfieldw = (1 - weight)*fieldw + weight*newfieldw;
            fieldt = dctI(tryfieldw);
            getpsi(tryfieldw);
            prop = prop + 1;
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
%     [tpsi, psisol] = ode45(@OCfpsi, [0 T], psi0, [], K, Vf, x, fieldt);
%     psisol = psisol.';
%     for ti = 1:(Nt+1)
%         psi(:, ti) = RK4interp(psisol, tpsi, (ti - 1)*dt);
%     end
%     evmiut = evmiu(psi, x);
%     evmiuw = dct(evmiut);
%    psi = deval(psisol, 0:dt:T);
    conv = conv(:, 1:(niter+1));
    niter
    weight
    prop
    J1 = sum(0.5*evmiuw.^2.*vfiltermiu)*dw;
    
    %%% Nested function: %%%
    function getpsi(fieldwint)
        [tpsi, psisol] = ode45(@OCfMLSpsi, [0 T], psi0, [], H0, miumat, fieldt, Nt, dt);
        psisol = psisol.';
        for ti = 1:(Nt+1)
            psi(:, ti) = RK4interp(psisol, tpsi, (ti - 1)*dt);
        end        
        evmiut = evmiuMLS(psi, miumat);
        evmiuw = dctI(evmiut);
        evmiufil = dctI(evmiuw.*vfiltermiu);
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J2fun = fieldwint(notsmall == 2).^2./vfilterE(notsmall == 2);
%        J2fun = J2fun(isfinite(J2fun));
 %       J2fun = J2fun(abs(J2fun)<100);
%        conv(:, niter) = (sum(0.5*evmiuw.^2.*vfiltermiu) - sum(fieldw(vfilterE~=0).^2./vfilterE(vfilterE~=0)))*dw;
        conv(:, niter + 1) = (sum(0.5*evmiuw.^2.*vfiltermiu) - sum(J2fun))*dw;
    end
end