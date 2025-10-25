function [fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, mniterc, J, J1, J2, Jorth, Jpnorm] = guessresults_pnaE0b(psi0, Vf, m,...
    xdomain, miux, a0, penalnormf, fguess, filterE, filtera, orthpenal, T, dt, Nt_ts, Nkr, tol)
% The program finds an optimal field for frequency control. It's intended
% for a Hamiltonian that depends on the spatial variable x.
% It uses a relaxation process.
% This version includes a penalty term on the loss of the norm of psi(T).
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
    a0M = spdiags(a0, 0, Nx, Nx);
    tolprop = 1e-3*tol;
    allapsi = zeros(Nx, allt_lasti);
    evcomT = 0;
    psiT = zeros(Nx, 1);
    compsiT = zeros(Nx, 1);
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
    allevat = zeros(1, allt_lasti);
    evaw = zeros(1, Nt + 1);
    J = 0;
    J1 = 0;
    J2 = 0;
    Jorth = 0;
    Jpnorm = 0;
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    vfilterE = zeros(1, Nt + 1);
%    vfiltermiu = zeros(1, Nwmiu + 1);
    vfiltera = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltera(wi) = filtera((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilterE(1), vfilterE(2:Nt), 0.5*vfilterE(Nt + 1)])*dw;
    coswT = ones(1, Nt + 1);
    coswT(2:2:(Nt + 1)) = -1;
    dctfilterET = sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1), vfilterE(2:Nt).*coswT(2:Nt), 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw;    
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    allpsi = zeros(Nx, allt_lasti);
%    chiihterm = zeros(Nx, allt_lasti);
    integw = chebgridw(Nt, Nt_ts, dt);
    igweights = integw(Nt_ts:2*(Nt_ts - 1))/dt;
    % ig is internal grid (within the time step).
    normpsiT = 0;
    if sqrt(2/pi)*abs(sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)]))*dw > tol ||...
            sqrt(2/pi)*abs(sum([0.5*fieldw(1)*coswT(1), fieldw(2:Nt).*coswT(2:Nt), 0.5*fieldw(Nt + 1)*coswT(Nt + 1)]))*dw > tol
        fieldw = getfieldwcon(fieldw);
    end
    getpsi(fieldw);
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = dctI(fieldw)/dctfactor;
    evat = allevat(1:(Nt_ts - 1):allt_lasti);
%    evaw = dctIfrom_ig(allevat, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    allevmiut = real(evmiu(allpsi, miux));
    evmiut = allevmiut(1:(Nt_ts - 1):allt_lasti);
    evmiuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
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
%        allapsi = a0M*allpsi + allpsi*spdiags(allfield.', 0, allt_lasti, allt_lasti);
        allapsi = a0M*allpsi;
        % The time-dependent part of the acceleration operator is omitted,
        % since it contributes low frequencies only. This just increases the
        % background of the signal, and does not contain any interesting
        % information.
        for allti = 1:allt_lasti
% The values of the expectation value are supposed to be real:
            allevat(allti) = real(allpsi(:, allti)'*allapsi(:, allti));
        end
        evaw = dctIfrom_ig(allevat, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        notsmall = (abs(fieldwint)>10*eps) + (abs(vfilterE)>10*eps);
        J1fun = 0.5*evaw.^2.*vfiltera*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J2fun = -fieldwint.^2./vfilterE*dw;
        J2fun([1, Nt + 1]) = J2fun([1, Nt + 1])/2;
        J2fun = J2fun(notsmall == 2);
        psiT = allpsi(:, allt_lasti);
        compsiT = ifft(K.*fft(miux.*psiT)) - miux.*ifft(K.*fft(psiT));
        evcomT = psiT'*compsiT;
        normpsiT = psiT'*psiT;
        J1 = sum(J1fun);
        J2 = sum(J2fun);
        Jorth = 0.5*orthpenal*evcomT^2;
        Jpnorm = penalnormf(normpsiT);
        J = J1 + J2 + Jorth + Jpnorm;
    end

end