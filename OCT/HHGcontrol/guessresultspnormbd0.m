function [fieldt, fieldw, psi, evmiut, evmiuw, evmiutcut, evmiuwcut, mniterc, J, J1, J2, Jpnorm] = guessresultspnormbd0(psi0, Vf, m,...
    xdomain, miux, penalnormf, fguess, filterE, filtermiu, tcutfun, Dtcutfun, T, dt, Nt_ts, Nkr, tol)
% The program finds the guess result for frequency control. It's intended
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
% Dpenalnormf: a function handle of the form @(x). Represents the
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
    miuxM = spdiags(miux, 0, Nx, Nx);
    %maxNiter = 1e3;
    tolprop = 1e-3*tol;
    allmiupsi = zeros(Nx, allt_lasti);
%    miualpsi = zeros(Nx, Nt + 1);
    psiT = zeros(Nx, 1);
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
    evmiuallt = zeros(1, allt_lasti);
    evmiuallt = zeros(1, allt_lasti);
    evmiualltcut = zeros(1, allt_lasti);
    evmiuwcut = zeros(1, Nt + 1);
    J = 0;
    J1 = 0;
    J2 = 0;
    Jpnorm = 0;
    t = 0:dt:T;
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt = [kron(t(1:Nt), ones(1, Nt_ts - 1)), T] + [kron(ones(1, Nt), t_ts(1:(Nt_ts - 1))), 0];
    tcutallt = tcutfun(allt);
    Dtcutallt = Dtcutfun(allt);
    alligweights = chebweights_gen(Nt_ts, [0 dt], t_ts(2:Nt_ts));
    igweights = [alligweights(1, Nt_ts - 1)*2, alligweights(2:(Nt_ts - 1), Nt_ts - 1).']/dt;
    vfilterE = zeros(1, Nt + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Nx, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    igweights = integw(Nt_ts:2*(Nt_ts - 1))/dt;
    % ig is internal grid (within the time step).
    normpsiT = 0;
    getpsi(fieldw);
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = dctI(fieldw)/dctfactor;
    evmiuallt = real(evmiu(allpsi, miux));
    evmiutcut = evmiualltcut(1:(Nt_ts - 1):allt_lasti);
    evmiut = evmiuallt(1:(Nt_ts - 1):allt_lasti);
    evmiuw = dctIfrom_ig(evmiuallt, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    
    %%% Nested functions: %%%
    function getpsi(fieldwint)
        allfield = dctIintgrid(fieldwint, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        [allpsi, field1, mniterc] = solveOCkr(@ihfieldmiux, K, Vf, x, psi0, [0 T], Nt, Nt_ts, Nkr, tolprop, allfield, miux);
        summniterc = summniterc + mniterc;
        allmiupsi = miuxM*allpsi;
        for allti = 1:allt_lasti
% The values of the expectation value are supposed to be real:
            evmiuallt(allti) = real(allpsi(:, allti)'*allmiupsi(:, allti));
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
        psiT = allpsi(:, allt_lasti);
        normpsiT = psiT'*psiT;
        J1 = sum(J1fun);
        J2 = sum(J2fun);
        Jpnorm = penalnormf(normpsiT);
        J = J1 + J2 + Jpnorm;
    end


end