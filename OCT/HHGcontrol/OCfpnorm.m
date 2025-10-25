function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfpnorm(psi0, Vf, m,...
    xdomain, miux, penalnormf, Dpenalnormf, fguess, filterE, filtermiu, orthpenal, iweight, T, dt, Nt_ts, Nkr, tol, maxNiter)
% The program finds an optimal field for frequency control. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% It uses a relaxation process.
% This version includes a penalty term on the loss of the norm of psi(T).
% It is possible to use any functional form for miu(x).
% The program uses the Arnoldi approach in the propagator. Hence, it is
% suitable for a Hamiltonian with complex eigenvalues (and for absorbing
% boundaries).
% The chiT is determined in a way that (d<miu>(T)/dt)^2 will be minimal,
% to prevent ringing.
% psi0:  The initial state
% Vf: the potential energy of H0, a function handle of the form: @(x). It
% is possible to insert the vector of the potential itself instead.
% m: The mass
% xdomain: The domain of the equally spaced grid of x: [min_x max_x)
% miux: a vector that represents miu.
% penalnormf: a function handle of the form @(x). Represents the penalty on
% the norm of psi(T).
% Dpenalnormf: a function handle of the form @(x). Represents the
% derivative of the function penalnormf.
% fguess, filterE, filtermiu, are function handles of w (the angular freqency), 
% or row vectors, where the w values are: 0:pi/T:pi/dt
% fguess: the initial guess for the field.
% filterE: the function that filters the field in the frequncy domain,
% multiplied by 1/(penalty factor)
% filtermiu: the function that filters the dipole eigenvalues in the
% frequncy domain.
% tol: the tolerance of the optimization process.
    Nt = T/dt;
    dw = pi/T;
    dctfactor = T/(sqrt(Nt*pi));
    Nx = length(psi0);
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/(2*m);
    if length(Vf) == 1
        % If Vf is a function handle:
        conjVf = @(x) conj(Vf(x));
    else
        % If Vf is a vector:
        conjVf = conj(Vf);
    end
%    miux = miufunx(x);
    miuxM = spdiags(miux, 0, Nx, Nx);
    %maxNiter = 1e3;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    allmiupsi = zeros(Nx, allt_lasti);
%    miualpsi = zeros(Nx, Nt + 1);
    evcomT = 0;
    psiT = zeros(Nx, 1);
    compsiT = zeros(Nx, 1);
    summniterc = 0;
    allfield = zeros(1, allt_lasti);
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
    allevmiut = zeros(1, allt_lasti);
    evmiuw = zeros(1, Nt + 1);
    evmiufil = zeros(1, Nt + 1);
    J1fun = zeros(1, Nt + 1);
    chimiupsi = zeros(1, allt_lasti);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    vfilterE = zeros(1, Nt + 1);
%    vfiltermiu = zeros(1, Nwmiu + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Nx, allt_lasti);
    chiihterm = zeros(Nx, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    igweights = integw(Nt_ts:2*(Nt_ts - 1))/dt;
    % ig is internal grid (within the time step).
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = iweight;
    niter = 0;
    normpsiT = 0;
    figure
    stop = false;
    uicontrol('style', 'push', 'string', 'Stop', 'callback', @buttonStop, 'position', [0, 0, 60, 20]);
    uicontrol('style', 'push', 'string', 'Pause', 'callback', @EnterPrompt, 'position', [70, 0, 60, 20]);
    getpsi(fieldw);
    plot(0:niter, real(conv(1:(niter + 1))));
    nprop = 1;
    while relE>tol && niter<maxNiter && ~stop
%        tic
        get_chiihterm(allmiupsi);
        chiT = orthpenal*evcomT*compsiT + Dpenalnormf(normpsiT)*psiT;
        [allchi, field1, mniterc] = solveOCkr(@ihalltchimiux, K, conjVf, x, chiT, [T 0], Nt, Nt_ts, Nkr, tolprop,...
            allfield(allt_lasti:-1:1), miux, chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        for allti = 1:allt_lasti
            chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*(miux.*allpsi(:, allti)));
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
        %drawnow
        plot(0:niter, real(conv(1:(niter + 1))))
        drawnow
    end
    if niter==maxNiter
        display('The program has failed to achieve the desired tolerance.')
    elseif stop
        display('The optimization process was stopped by the user.')
    end
    conv = conv(1:(niter+1));
    niter
    weight
    nprop
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = dctI(fieldw)/dctfactor;
%    allevmiut = evmiu(allpsi, miux);
    allevmiut = real(evmiu(allpsi, miux));
    evmiut = allevmiut(1:(Nt_ts - 1):allt_lasti);
    evmiuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    mallniterc = summniterc/nprop;
    J1 = sum(J1fun);
    for allti = 1:allt_lasti
        chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*(miux.*allpsi(:, allti)));
    end
    grad = zeros(1, Nt + 1);
    grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
    grad = grad + 2*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    maxgrad = max(abs(grad));
    evorth = 1i*evcomT
    beep
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCkr(@ihfieldmiux, K, Vf, x, psi0, [0 T], Nt, Nt_ts, Nkr, tolprop, allfield, miux);
        summniterc = summniterc + mniterc;
        allmiupsi = miuxM*allpsi;
        for allti = 1:allt_lasti
%            allevmiut(allti) = allpsi(:, allti)'*allmiupsi(:, allti);
% The values of the expectation value are supposed to be real:
            allevmiut(allti) = real(allpsi(:, allti)'*allmiupsi(:, allti));
        end
        evmiuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evmiuw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        psiT = allpsi(:, allt_lasti);
        compsiT = ifft(K.*fft(miux.*psiT)) - miux.*ifft(K.*fft(psiT));
        evcomT = psiT'*compsiT;
        normpsiT = psiT'*psiT;
        conv(niter + 1) = sum(J1fun) + sum(J2fun) + 0.5*orthpenal*evcomT^2 + penalnormf(normpsiT);
    end

    function get_chiihterm(allmiupsi)
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        chiihterm = -allmiupsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1).', 0, allt_lasti, allt_lasti);     
    end

    function buttonStop(hObject, event)
        stop = true;
    end

end

function EnterPrompt(hObject, event)
    keyboard
end