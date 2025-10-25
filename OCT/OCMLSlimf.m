function [fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimf(psi0, target, E0, Edomain, miu, fguess, ...
    filterE, T, dt, Nt_ts, Ncheb, tol)
% The program finds an optimal field for a state to state transition, 
% with a limitation on the allowed frequency.
% It's intended for a many level system.
% It uses a relaxation process.
% psi0:  The initial state
% target: The target state
% E0: The eigenenergies of H0
% miu: The dipole operator matrix
% fguess, filterE, are function handles of w (the angular freqency), 
% or row vectors, when the w values are: 0:pi/T:pi/dt
% fguess: the initial guess for the field.
% filterE: the function that filters the field in the frequncy domain,
% multiplied by 1/(penalty factor)
% tol: the tolerance of the optimization process.
    Nt = T/dt;
    dw = pi/T;
    dctfactor = T/(sqrt(Nt*pi));
%    Nwfield = ceil(maxwE/dw);
%    Nwmiu = ceil(maxwmiu/dw);
    Npsi = length(psi0);
    maxNiter = 300;
    tolprop = 1e-3*tol;
    conv = zeros(1, maxNiter + 1);
    H0 = diag(E0);
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
    J1 = 0;
    chimiupsi = zeros(1, Nt + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    psi = zeros(Npsi, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Npsi, Nt_ts, Nt);
    overlap = 0;
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = 1;
    niter = 0;
    getpsi(fieldw);
    nprop = 1;
    while relE>tol && niter<maxNiter
%        tic
        chiT = overlap*target;
        [allchi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tolprop,...
            allfield(allt_lasti:-1:1));
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
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    niter
    weight
    nprop
    fieldt = dctI(fieldw)/dctfactor;
    mallniterc = summniterc/nprop;
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
        psiT = allpsi(:, Nt_ts, Nt);
        overlap = target'*psiT;
        J1 = overlap*conj(overlap);
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        conv(niter + 1) = J1 + sum(J2fun);
    end

end