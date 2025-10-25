function [fieldt, fieldw, psi, relE, conv, niter, mallniterc, J1, maxgrad] = OCMLSlimfE0b(psi0, target, E0, Edomain, miu, fguess, ...
    filterE, T, dt, Nt_ts, Ncheb, tol)
% The program finds an optimal field for a state to state transition, 
% with a limitation on the allowed frequency.
% It's intended for a many level system.
% It uses a relaxation process.
% The computation of fieldw is more accurate than in OCMLSlimf.
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
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    chimiupsi = zeros(1, allt_lasti);
    igweights = chebweights(Nt_ts, 1);
    igweights = [igweights(1)*2, igweights(2:(Nt_ts - 1))];
    psi = zeros(Npsi, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    allpsi = zeros(Npsi, Nt_ts, Nt);
    dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilterE(1), vfilterE(2:Nt), 0.5*vfilterE(Nt + 1)])*dw;
    coswT = ones(1, Nt + 1);
    coswT(2:2:(Nt + 1)) = -1;
    dctfilterET = sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1), vfilterE(2:Nt).*coswT(2:Nt), 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw;    
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    overlap = 0;
    lastfieldw = fieldw;
    relE = tol + 1;
    weight = 1;
    niter = 0;
    if sum([0.5*fieldw(1), fieldw(2:Nt), 0.5*fieldw(Nt + 1)]) > tol ||...
            sum([0.5*fieldw(1)*coswT(1), fieldw(2:Nt).*coswT(2:Nt), 0.5*fieldw(Nt + 1)*coswT(Nt + 1)]) > tol
        fieldw = getfieldwcon(fieldw);
    end
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
            for ti = 1:(Nt_ts - 1)
                chimiupsi((tsi - 1)*(Nt_ts - 1) + ti) = -imag(allchi(:, Nt_ts - ti + 1, Nt - tsi + 1)'*miu*allpsi(:, ti, tsi));
            end
        end
        chimiupsi(allt_lasti) = -imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt));
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
        for ti = 1:(Nt_ts - 1)
            chimiupsi((tsi - 1)*(Nt_ts - 1) + ti) = -imag(allchi(:, Nt_ts - ti + 1, Nt - tsi + 1)'*miu*allpsi(:, ti, tsi));
        end
    end
    chimiupsi(allt_lasti) = -imag(allchi(:, 1, 1)'*miu*allpsi(:, Nt_ts, Nt));
    grad = zeros(1, Nt + 1);
    grad(fieldw ~= 0) = 2*fieldw(fieldw ~= 0)./vfilterE(fieldw ~= 0);
    grad = grad + 2*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    maxgrad = max(abs(grad));
    
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