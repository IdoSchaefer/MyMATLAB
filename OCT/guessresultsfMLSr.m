function [fieldt, fieldw, psi, evmiut, evmiuw, mniterc, J, J1, J2, Jorth] = guessresultsfMLSr(psi0, E0, Edomain, miu, ...
    fguess, filterE, filtermiu, orthpenal, T, dt, Nt_ts, Ncheb, tol)
% The program finds the guess results for frequency control. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% The chiT is determined in a way that (d<miu>(T)/dt)^2 will be minimal,
% to prevent ringing.
% psi0:  The initial state
% Vf: the potential energy of H0, a function handle of the form: @(x)
% m: The mass
% xdomain: The domain of the equally spaced grid of x: [min_x max_x)
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
    Npsi = length(psi0);
    tolprop = 1e-3*tol;
    H0 = diag(E0);
    % comH0miu = [H, miu]:
    comH0miu = H0*miu - miu*H0;
    evcomT = 0;
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
    J1 = 0;
    J2 = 0;
    Jorth = 0;
%     chimiupsi = zeros(1, Nt + 1);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
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
    getpsi(fieldw);
    fieldt = dctI(fieldw)/dctfactor;
%     for tsi = 1:Nt
%         chimiupsi(tsi) = -imag(allchi(:, Nt_ts, Nt - tsi + 1)'*(x.*allpsi(:, 1, tsi)));
%     end
%     chimiupsi(Nt + 1) = -imag(allchi(:, 1, 1)'*(x.*allpsi(:, Nt_ts, Nt)));
%     grad = zeros(1, Nt + 1);
%     grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
%     grad = grad + 2*dctI(chimiupsi)*dctfactor;
%     maxgrad = max(abs(grad));
    evorth = 1i*evcomT
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        tic
       [allpsi, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tolprop, allfield);
        toc
        psi(:, 1:Nt) = allpsi(:, 1, :);
        psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
        evmiut = evmiuMLS(psi, miu);
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
        evcomT = psiT'*comH0miu*psiT;
        J1 = sum(J1fun);
        J2 = sum(J2fun);
        Jorth = 0.5*orthpenal*evcomT^2;
        J = J1 + J2 + Jorth;
    end

end