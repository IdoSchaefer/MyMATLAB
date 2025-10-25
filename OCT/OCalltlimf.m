function [fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCalltlimf(psi0, target, Vf, m, Edomain,...
    xdomain, coeft, coefT, fguess, filterE, iweight, T, dt, Nt_ts, Ncheb, tol)
% The program is intended to find an optimal field for a state to state transition, as fast as possible.
% For this purpose, there is the following term in the functional:
% coeft*integ(<psi(t)|target><target|psi(t)>), in addition to the term in
% the final time, multiplied by coefT.
% The program includes a limitation on the allowed frequency. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% It uses a relaxation process.
% psi0:  The initial state
% target: The target state
% Vf: the potential energy of H0, a function handle of the form: @(x)
% xdomain: The domain of the equally spaced grid of x: [min_x max_x)
% fguess, filterE, are function handles of w (the angular freqency), 
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
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/(2*m);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    summniterc = 0;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
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
    J1 = 0;
    chimiupsi = zeros(1, Nt + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    psi = zeros(Nx, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Nx, Nt_ts, Nt);
    overlaps = zeros(1, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    %chiihterm = zeros(Nx, Nt_ts, Nt);
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = iweight;
    niter = 0;
    figure
    getpsi(fieldw);
    plot(0:niter, real(conv(1:(niter + 1))));
    nprop = 1;
    while relE>tol && niter<maxNiter
%        tic
%        get_chiihterm();
        chiihterm = -coeft*target*overlaps(allt_lasti:-1:1);
        chiT = -chiihterm(:, 1)/coeft*coefT;
        [allchi, field1, mniterc] = solveOCih(@ihalltchi, K, Vf, Edomain, x, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1), chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 1;
        for tsi = 1:Nt
            chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*(x.*allpsi(:, 1, tsi)));
        end
        chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*(x.*allpsi(:, Nt_ts, Nt)));
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
        plot(0:niter, real(conv(1:(niter + 1))))
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    conv = conv(1:(niter+1));
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    niter
    weight
    nprop
    fieldt = dctI(fieldw)/dctfactor;
    mallniterc = summniterc/nprop;
    for tsi = 1:Nt
        chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*(x.*allpsi(:, 1, tsi)));
    end
    chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*(x.*allpsi(:, Nt_ts, Nt)));
    grad = zeros(1, Nt + 1);
    grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
    grad = grad + 2*dctI(chimiupsi)*dctfactor;
    maxgrad = max(abs(grad));
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCih(@ihfield, K, Vf, Edomain, x, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        summniterc = summniterc + mniterc;
        for ti = 1:(Nt_ts - 1)
            psiti(:, :) = allpsi(:, ti, :);
            overlaps(ti:(Nt_ts - 1):(allt_lasti - 1)) = target'*psiti;
        end
        overlaps(allt_lasti) = target'*allpsi(:, Nt_ts, Nt);
        J1T = coefT*overlaps(allt_lasti)*conj(overlaps(allt_lasti));
        J1t = coeft*sum(integw.*overlaps.*conj(overlaps));
        J1 = J1T + J1t;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        conv(niter + 1) = J1 + sum(J2fun);
    end

end