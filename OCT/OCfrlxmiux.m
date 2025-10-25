function [fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfrlxmiux(psi0, Vf, m, Edomain, xdomain,...
    miufunx, fguess, filterE, filtermiu, orthpenal, iweight, T, dt, Nt_ts, Ncheb, tol)
% The program finds an optimal field for frequency control. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% It uses a relaxation process.
% The chiT is determined in a way that (d<miu>(T)/dt)^2 will be minimal,
% to prevent ringing.
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
    miux = miufunx(x);
    miuxM = spdiags(miux, 0, Nx, Nx);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    allmiupsi = zeros(Nx, allt_lasti);
    miupsi = zeros(Nx, Nt + 1);
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
    evmiut = zeros(1, Nt + 1);
    evmiuw = zeros(1, Nt + 1);
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
%    chiihterm = zeros(Nx, Nt_ts, Nt);
    chiihterm = zeros(Nx, allt_lasti);
    lastfieldw = fieldw;
    relE = tol + 1;
%    weight = 1;
    weight = iweight;
    niter = 0;
    getpsi(fieldw);
    nprop = 1;
    while relE>tol && niter<maxNiter
%        tic
        get_chiihterm(allpsi);
        chiT = orthpenal*evcomT*compsiT;
        % The function ihalltchimiux uses another form of chiihterm. This program has to be corrected accordingly.
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
    evorth = 1i*evcomT
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCih(@ihfieldmiux, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield, miux);
        summniterc = summniterc + mniterc;
        psi(:, 1:Nt) = allpsi(:, 1, :);
        psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
        miupsi = miuxM*psi;
        for tsi = 1:(Nt + 1)
            evmiut(tsi) = psi(:, tsi)'*miupsi(:, tsi);
        end
%        evmiut = evmiu(psi, miux);
        evmiuw = dctI(evmiut)*dctfactor;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evmiuw.^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        %integrand = (-0.5*evmiuw.^2.*vfiltermiu,  + J2fun)*dw/(Nt_ts-1);
        %integrand([1, Nt*(Nt_ts - 1)]) = integrand([1, Nt*(Nt_ts-1)])/2;
        %minusJ = sum(integrand); 
        psiT = psi(:, Nt + 1);
        compsiT = ifft(K.*fft(miux.*psiT)) - miux.*ifft(K.*fft(psiT));
        evcomT = psiT'*compsiT;
        conv(niter + 1) = sum(J1fun) + sum(J2fun) + 0.5*orthpenal*evcomT^2;
    end

    function get_chiihterm(allpsi)
        evmiufil = dctIintgrid(evmiuw.*vfiltermiu, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        allmiupsi(:, 1:(Nt_ts - 1):allt_lasti) = miupsi;
        for ti = 2:(Nt_ts - 1)
            psiti(:, :) = allpsi(:, ti, :);
            allmiupsi(:, ti:(Nt_ts - 1):(allt_lasti - 1)) = miuxM*psiti;
        end
        allmiupsi(:, allt_lasti) = miuxM*allpsi(:, Nt_ts, Nt);
        chiihterm = -allmiupsi(:, allt_lasti:-1:1)*spdiags(evmiufil(allt_lasti:-1:1).', 0, allt_lasti, allt_lasti); 
        
%         for ti = 1:Nt_ts
%             chiihterm(:, ti, 1) = -evmiufil(Nt*(Nt_ts - 1) - ti + 2)*miux.*allpsi(:, Nt_ts - ti + 1, Nt);
%         end        
%         for tsi = 2:Nt
%             chiihterm(:, 1, tsi) = chiihterm(:, Nt_ts, tsi - 1);
%             for ti = 2:Nt_ts
%                 chiihterm(:, ti, tsi) = -evmiufil((Nt - tsi + 1)*(Nt_ts - 1) - ti + 2)*miux.*allpsi(:, Nt_ts - ti + 1, Nt - tsi + 1);
%             end
%         end    
      end

end