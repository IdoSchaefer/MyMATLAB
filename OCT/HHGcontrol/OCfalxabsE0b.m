function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfalxabsE0b(psi0, Vf, m,...
    xdomain, miufunx, allowedf, penalforbf, fguess, filterE, filtermiu, orthpenal, iweight, T, dt, Nt_ts, Nkr, tol)
% The program finds an optimal field for frequency control. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% It uses a relaxation process.
% This version includes the option of restricting the allowed x zone.
% It is possible to use any functional form for miu(x).
% The program uses the Arnoldi approach in the propagator. Hence, it is
% suitable for a Hamiltonian with complex eigenvalues (and for absorbing
% boundaries).
% The chiT is determined in a way that (d<miu>(T)/dt)^2 will be minimal,
% to prevent ringing.
% psi0:  The initial state
% Vf: the potential energy of H0, a function handle of the form: @(x)
% m: The mass
% xdomain: The domain of the equally spaced grid of x: [min_x max_x)
% miufunx: a function of x that represents miu.
% allowedf: a function of x that determines the contribution of psi(x)
% to J1.
% penalforbf: a function of x; the penalty of psi(x).
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
    if length(miufunx) == 1
        % If miufunx is a function handle:
        miux = miufunx(x);
    else
        % If miufunx is a vector:
        miux = miufunx;
    end
    miuxal = miux.*allowedf(x).^2;
    miuxalM = spdiags(miuxal, 0, Nx, Nx);
    penalforbv = penalforbf(x);
    projforb = spdiags(penalforbv, 0, Nx, Nx);
    maxNiter = 1e3;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    allmiualpsi = zeros(Nx, allt_lasti);
%    miualpsi = zeros(Nx, Nt + 1);
    allpfpsi = zeros(Nx, allt_lasti);
    allevprojforb = zeros(1, allt_lasti);
    evcomT = 0;
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
    allevmiualt = zeros(1, allt_lasti);
    evmiualw = zeros(1, Nt + 1);
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
    dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilterE(1), vfilterE(2:Nt), 0.5*vfilterE(Nt + 1)])*dw;
    coswT = ones(1, Nt + 1);
    coswT(2:2:(Nt + 1)) = -1;
    dctfilterET = sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1), vfilterE(2:Nt).*coswT(2:Nt), 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw;    
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    allpsi = zeros(Nx, allt_lasti);
    chiihterm = zeros(Nx, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    igweights = integw(Nt_ts:2*(Nt_ts - 1))/dt;
    % ig is internal grid (within the time step).
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = iweight;
    niter = 0;
    figure
    stop = false;
    uicontrol('style', 'push', 'string', 'Stop', 'callback', @buttonStop, 'position', [0, 0, 60, 20]);
    uicontrol('style', 'push', 'string', 'Pause', 'callback', @EnterPrompt, 'position', [70, 0, 60, 20]);
    if sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)]) > tol ||...
            sum([0.5*fieldw(1)*coswT(1), fieldw(2:Nt).*coswT(2:Nt), 0.5*fieldw(Nt + 1)*coswT(Nt + 1)]) > tol
        fieldw = getfieldwcon(fieldw);
    end
    getpsi(fieldw);
    plot(0:niter, real(conv(1:(niter + 1))));
    nprop = 1;
    while relE>tol && niter<maxNiter && ~stop
%        tic
        get_chiihterm(allmiualpsi);
        chiT = orthpenal*evcomT*compsiT;
        [allchi, field1, mniterc] = solveOCkr(@ihalltchimiux, K, conjVf, x, chiT, [T 0], Nt, Nt_ts, Nkr, tolprop,...
            allfield(allt_lasti:-1:1), miux, chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        for allti = 1:allt_lasti
            chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*(miux.*allpsi(:, allti)));
        end
        niter = niter + 1;
        flag = 0;
        newfieldw_unc = dctfactor*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights).*vfilterE;
        newfieldw = getfieldwcon(newfieldw_unc);
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
        plot(0:niter, real(conv(1:(niter + 1))))
        drawnow
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
    allevmiut = evmiu(allpsi, miux);
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
    
    function fieldw_con = getfieldwcon(fieldw_unc)
    % The function computes a constrained field spectrum, fieldw_con, with 0 boundaries in the
    % time domain, from the unconstrained field spectrum fieldw_unc.
        fieldt_unc0 = sqrt(2/pi)*sum([0.5*fieldw_unc(1), fieldw_unc(2:Nt), 0.5*fieldw_unc(Nt + 1)])*dw;    
        fieldt_uncT = sqrt(2/pi)*...
            sum([0.5*fieldw_unc(1)*coswT(1), fieldw_unc(2:Nt).*coswT(2:Nt), 0.5*fieldw_unc(Nt + 1)*coswT(Nt + 1)])*dw;    
        lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
        lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;
        fieldw_con = fieldw_unc - vfilterE.*(lambda0 + lambdaT*coswT);
    end
    
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCkr(@ihfieldmiux, K, Vf, x, psi0, [0 T], Nt, Nt_ts, Nkr, tolprop, allfield, miux);
        summniterc = summniterc + mniterc;
        allmiualpsi = miuxalM*allpsi;
        allpfpsi = projforb*allpsi;
        for allti = 1:allt_lasti
            allevmiualt(allti) = allpsi(:, allti)'*allmiualpsi(:, allti);
            allevprojforb(allti) = allpsi(:, allti)'*allpfpsi(:, allti);
        end
        evmiualw = dctIfrom_ig(allevmiualt, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evmiualw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
%         for ti = 1:(Nt_ts - 1)
%             psiti(:, :) = allpsi(:, ti, :);
%             allpfpsi(:, ti:(Nt_ts - 1):(allt_lasti - 1)) = projforb*psiti;
%             for tsi = 1:Nt
%                 allevprojforb((tsi - 1)*(Nt_ts - 1) + ti) = allpsi(:, ti, tsi)'*allpfpsi(:, (tsi - 1)*(Nt_ts - 1) + ti);
%             end
%         end
%         allpfpsi(:, allt_lasti) = projforb*allpsi(:, Nt_ts, Nt);
%         allevprojforb(allt_lasti) = allpsi(:, Nt_ts, Nt)'*allpfpsi(:, allt_lasti);
        psiT = allpsi(:, allt_lasti);
        compsiT = ifft(K.*fft(miuxal.*psiT)) - miuxal.*ifft(K.*fft(psiT));
        evcomT = psiT'*compsiT;
        conv(niter + 1) = sum(J1fun) + sum(J2fun) - sum(integw.*allevprojforb) + 0.5*orthpenal*evcomT^2;
    end

    function get_chiihterm(allmiualpsi)
        evmiufil = dctIintgrid(evmiualw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        chiihterm = allpfpsi(:, allt_lasti:-1:1) -...
            allmiualpsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1).', 0, allt_lasti, allt_lasti); 
        
      end

    function buttonStop(hObject, event)
        stop = true;
    end

end

function EnterPrompt(hObject, event)
    keyboard
end