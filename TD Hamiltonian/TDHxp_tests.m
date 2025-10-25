function [U, mniter, matvecs, errors, Ferror, texp_dif] = TDHxp_tests(K, V, Vtfun, ihfun, Edomain, ui, x, tdomain, Nts, Nt_ts, Ncheb, tol, Niter, Niter1st, varargin)
% The program solves Schrodinger equation for a time dependent Hamiltonian.
% It uses the Hillel Tal-Ezer's propagator.
% The Hamiltonian can be nonlinear (u dependent).
% There is an option to include an inhomogeneos term.
% The program should be used if H0 is composed only from x diagonal and p
% diagonal elements, and the time dependent perturbation Vt is composed
% only from x diagonal elements.
% Input:
% K: the p diagonal element of the unperturbed Hamiltonian (kinetic
% energy).
% V: the x diagonal element of the unperturbed Hamiltonian (potential
% energy).
% Vtfun: A function handle of the form: @(u, x, t, more arguments). Returns
% the x vector of the time dependent, nonlinear perturbation: Vt(u, x, t).
% ihfun: A function handle of the form: @(t, more arguments). Returns the
% inhomogeneous term vector. If there is no inhomogeneous term, put []
% instead.
% Note: "more arguments" must be the same for the 2 functions. The arguments
% will be written in the place of "varargin" in this program.
% Edomain: the domain of the eigenenergies of the hamiltonian: [minE maxE]
% ui: the initial state vector.
% x: the x grid
% tdomain: the time interval of the propagation: [ti tf]
% Nts: the number of time steps. 
% Nt_ts: the number of interior Chebyshev time points in the time step, used during
% the computational process.
% Ncheb: the number of terms used to expand the function of matrix in the
% Chebyshev expansion.
% tol: the desired tolerance of the convergence.
% Output:
% U: contains the solution in the points: ti:((tf - ti)/Nts):tf, in different columns.
% mniter: the mean number of iteration for a time step; should be close to
% 1 for ideal efficiency of the algorithm.
% matvecs: the number of Hamiltonian operations.
    tinit = tdomain(1);
    tf = tdomain(2);
    % The length of the time interval of the whole propagation(can be negative):
    T = tf - tinit;
    % The length of the time step interval:
    Tts = T/Nts;
    % The index of the middle term in the time step:
    tmidi = round(Nt_ts/2);
    minE = Edomain(1);
    maxE = Edomain(2);
    Nx = length(ui);
    U = zeros(Nx, Nts + 1);
    U(:, 1) = ui;
    % The Chebyshev points for expansion in time, in the domain in which the Chebyshev expansion
    % is defined: [-1 1]
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The Chebyshev points for expansion in time, in the domain of the time
    % variable (changing of variables):
    t_ts = 0.5*(tcheb+1)*Tts;
    test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
%    test_tpoint = (t_ts(Nt_ts - 1) + t_ts(Nt_ts))/2;
    tp_theta = acos((tcheb(tmidi) + tcheb(tmidi + 1))/2);
    texpEfactor = -2*sin((Nt_ts - 1)*tp_theta)*(cos(2*tp_theta) - 1)/sin(tp_theta);
    abs_texpEfactor = abs(texpEfactor);
    texpEfactor1 = 8*Tts/((Nt_ts - 1)*((Nt_ts - 1)^2 - 4));
    texpEfactor2 = abs(4*Tts^2*(Nt_ts - 1)/((Nt_ts - 1)^4 - 10*(Nt_ts - 1)^2 + 9));
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation:
%    t_2ts = [t_ts.'; Tts + t_ts(2:Nt_ts).'].';
    t_2ts = [t_ts, test_tpoint, Tts + t_ts(2:Nt_ts), Tts + test_tpoint];
    Ninteg = 60;
    tts_integ = Tts/(2*(Ninteg)):Tts/(Ninteg):(Tts - Tts/(2*(Ninteg)));
    tint_2ts = [tts_integ, Tts + tts_integ];
    tcheb_int = 2*tts_integ/Tts - 1;
    theta_int = acos(tcheb_int);
    % Computing the coefficients for the Chebyshev expansion of the
    % function of matrix, in all the interior time points.
    % CchebFts contains the coefficients of the current time step, and
    % CchebFnext contains the coefficients of the next one.
%    [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb);
%    [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb, tol);
    [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb, eps);
    CchebF_int = FchebC1(tint_2ts, Nt_ts, minE, maxE, Ncheb, eps);
    CchebFts_int = CchebF_int(:, 1:Ninteg);
    CchebFnext_int = CchebF_int(:, (Ninteg + 1):(2*Ninteg));
    EdomainL = maxE - minE;
    xtest = ((minE + EdomainL/(Nx + 1)):EdomainL/(Nx + 1):(maxE - EdomainL/(Nx + 1))).';
%    Fxtest = Ffun(xtest, Tts, Nt_ts, Nx);
%    Fxtest = Ffun(xtest, Tts, Nt_ts, tol);
    Fxtest = Ffun(xtest, Tts, Nt_ts, eps);
    Ferror = max(abs(chebc2result(CchebFts(:, Nt_ts), Edomain, xtest) - Fxtest)./abs(Fxtest));
    % Computing the matrix of the Taylor polynomials in time.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    [timeMts, timeMnext] = maketM(t_2ts, Nt_ts);
    timeM_int = maketM1(tint_2ts, Nt_ts);
    timeMts_int = timeM_int(:, 1:Ninteg);
    timeMnext_int = timeM_int(:, (Ninteg + 1):(2*Ninteg));    
    % Computing the coefficients of the transformation from the Newton
    % interpolation polynomial terms, to a Taylor like form:
    Cr2t = r2Taylor4(t_ts, Tts);
    % The solution in the interior points in the time step.
    % Every column represents an interior time point in the time step:
    %Unew = zeros(Nx, Nt_ts);
    % The "inhomogeneos" vectors, including the term that is realy u
    % dependent (appears in all the cases, not only for the inhomogeneous Schrodinger equation):
%    Fi = zeros(Nx, Nt_ts);
    Fi = zeros(Nx, Nt_ts + 1);
    Fi_new = zeros(Nx, Nt_ts);
    Fi_int = zeros(Nx, Nt_ts);
    Cnewton = zeros(Nx, Nt_ts + 1);
    Cnewton_der = zeros(Nx, 2);
    % The Lambda vectors are defined in a recursive way, and contain information about
    % the time dependence of the Fi vectors:
    Lambda = zeros(Nx, Nt_ts+1);
    % If there is no inhomogeneous term in the equation, ihfun = [], and the length is 0 (false).
    % If there is, the length is 1 (true).
    there_is_ih = length(ihfun);
    if there_is_ih
%        ihterm = zeros(Nx, Nt_ts);
        ihterm = zeros(Nx, Nt_ts + 1);
        ihterm(:, 1) = ihfun(0, varargin{:});
        ihterm_integ = zeros(Nx, Ninteg);
    end
    DM = DchebbMat(Nt_ts - 1, Tts);
    DM = DM(Nt_ts:-1:1, Nt_ts:-1:1);
    % The 0'th order approximation is the first guess, for the first time step.
    % Every column represents an interior time point in the time step:
%    Uguess = guess_ts1(ui, Nt_ts);
    errors = zeros(Nts, 18);
    texp_dif = zeros(Nx, Nts);
    Uguess = guess_ts1(ui, Nt_ts + 1);
    Uguess_int = guess_ts1(ui, Ninteg);
    allniter = 0;
    if nargin<13
        Niter = 10;
        Niter1st = 13;
    end
    Cnewton_int = zeros(Nx, Ninteg);
    texpansion_error7 = NaN;
%     texpansion_error8 = NaN;
    wcheb = chebweights(Nt_ts, Tts).';
    for tsi = 1:Nts
        % The time, represented by the interior time points:
%        t = tinit + Tts*(tsi - 1) + t_ts;
        t = tinit + Tts*(tsi - 1) + [t_ts, test_tpoint];
        t_integ = tinit + Tts*(tsi - 1) + tts_integ;
        % The first guess for the iterative process, for the convergence of the u
        % values. Every column represents an interior time point in the time step:
        Ulast = Uguess;
        Ulast_int = Uguess_int;
        Lambda(:, 1) = U(:, tsi);
        if there_is_ih
            % Computing the inhomogeneous term:
%            for ti = 2:Nt_ts
            for ti = 2:(Nt_ts + 1)
                ihterm(:, ti) = ihfun(t(ti), varargin{:});
            end
            for ti = 1:Ninteg
                ihterm_integ(:, ti) = ihfun(t_integ(ti), varargin{:});
            end            
        end
        % Starting an iterative process, until convergence:
        niter = 0;
        reldif = tol + 1;
%         Lambda(:, Nt_ts + 1) = ui;
%         Vts = V;
        while (reldif>tol && ((tsi>1 && niter<Niter) || (tsi == 1 && niter<Niter1st)))
            % Vthalf is the x vector representing Vt(t = t(tmidi)): 
            Vthalf = Vtfun(Ulast(:, tmidi), x, t(tmidi), varargin{:});
            % Vts is the x vector representing the total potential energy 
            % operator at t = t(tmidi):
            %[Vkr, Hkr] = createKrop(@(u) Hpsi(K, Vts, u), Lambda(:, Nt_ts + 1), Ncheb);
            Vts = V + Vthalf;
            % Calculation of the inhomogeneous fi vectors:
%            for ti = 1:Nt_ts
            for ti = 1:(Nt_ts + 1)
                Fi(:, ti) = -1i*(Vtfun(Ulast(:, ti), x, t(ti), varargin{:}) - Vthalf).*Ulast(:, ti);
            end
            for ti = 1:(Ninteg)
                Fi_int(:, ti) = -1i*(Vtfun(Ulast_int(:, ti), x, t_integ(ti), varargin{:}) - Vthalf).*Ulast_int(:, ti);
            end
            % If there is an inhomogeneous term, we add it to the fi
            % vectors:
            if there_is_ih
                Fi = Fi + ihterm;
                Fi_int = Fi_int + ihterm_integ;
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation in the points t_ts.
            % The divided differences are computed by the function divdif.
            % For numerical stability, we have to transform the time points
            % in the time step, to points in an interval of length 4:
%             [Cnewton, diagonal] = divdif([t_ts, test_tpoint]*4/Tts, Fi);
            [Cnewton(:, 1:Nt_ts), diagonal] = divdif(t_ts*4/Tts, Fi(:, 1:Nt_ts));
            Cnewton(:, Nt_ts + 1) = new_divdif([t_ts, test_tpoint]*4/Tts, Fi(:, Nt_ts + 1), diagonal);
            for tinti = 1:Ninteg
%                 Cnewton_temp = divdif([t_ts, tts_integ(tinti)]*4/Tts, [Fi(:, 1:Nt_ts), Fi_int(:, tinti)]);
%                 Cnewton_int(:, tinti) = Cnewton_temp(:, Nt_ts + 1);
                Cnewton_int(:, tinti) = new_divdif([t_ts, tts_integ(tinti)]*4/Tts, Fi_int(:, tinti), diagonal);
            end
            for tderi = 1:2
%                 Cnewton_temp = divdif([t_ts, tts_integ(tderi - 1 + Ninteg/2)]*4/Tts, [Fi(:, 1:Nt_ts), Fi_int(:, tderi - 1 + Ninteg/2)]);
%                 Cnewton_der(:, tderi) = Cnewton_temp(:, Nt_ts + 1);
                Cnewton_der(:, tderi) = new_divdif([t_ts, tts_integ(tderi - 1 + Ninteg/2)]*4/Tts, Fi_int(:, tderi - 1 + Ninteg/2), diagonal);
            end
            % Calculating the Taylor like coefficients:
            Ctaylor = zeros(Nx, Nt_ts);
            Ctaylor(:, 1) = Cnewton(:, 1);
            for Newtoni = 2:Nt_ts
                Ctaylor(:, 1:Newtoni) = Ctaylor(:, 1:Newtoni) + Cnewton(:, Newtoni)*Cr2t(Newtoni, 1:Newtoni);
            end
            % Calculation of the Lambda vectors:
            for polyi = 2:(Nt_ts+1)
                Lambda(:, polyi) = (-1i*Hpsi(K, Vts ,Lambda(:, polyi-1)) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            % Vcheb is a matrix. Contains the vectors:
            % Tn(Hts)*Lambda(: ,Nt_ts + 1),  n = 0, 1, ..., Ncheb-1
            % where the Tn(x) are the Chebyshev polynomials.
            % The n'th vector is the (n+1)'th column of Vcheb.
            % The vectors will be used to approximate the function of
            % the operator Hts = K + Vts, that operates on Lambda(:, Nt_ts+1).
            Vcheb = vchebM(K, Vts, Lambda(:, Nt_ts + 1), Nx, minE, maxE, Ncheb);
            % Calculation of the wave function in all the time points
            % within the time step:
            [Unew, nFnU] = UfromLamb(Lambda, timeMts, Vcheb, CchebFts);
            Unew_int = UfromLamb(Lambda, timeMts_int, Vcheb, CchebFts_int);
            % To check the convergence of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            DUnew = Unew(:, 1:Nt_ts)*DM.';
            reldifD = norm(DUnew(:, tmidi) + 1i*Hpsi(K, Vts ,Unew(:, tmidi)))/norm(DUnew(:, tmidi));
%            texpansion_error = 2*norm(Cnewton(:, Nt_ts + 1))/norm(Fi(:, Nt_ts + 1));
            texpansion_error = abs_texpEfactor*norm(Cnewton(:, Nt_ts + 1))/norm(Fi(:, Nt_ts + 1));
            texpE_RHS = abs_texpEfactor*norm(Cnewton(:, Nt_ts + 1))/norm(-1i*Hpsi(K, Vts, Unew(:, Nt_ts + 1)) +  Fi(:, Nt_ts + 1));
            texpansion_error10 = abs_texpEfactor*norm(Cnewton(:, Nt_ts + 1))/norm(-1i*Hpsi(K, Vts, Unew(:, tmidi)) +  Fi(:, tmidi))*norm(Unew(:, Nt_ts) - Unew(:, 1))/norm(Unew(:, 1));
            % tmidi can be safely replaced by 1. Note that -1i*Hpsi(K, Vts, Unew(:, 1)) is already computed in the computation of 
            % Lambda(:, 2)
            % For a nonlinear problem, Vthalf should be fixed.
            % The following estimations yield identical results to
            % texpansion_error and texpE_RHS:
%             texpansion_error2 = norm(Fi_testpoint - Fi(:, Nt_ts + 1))/norm(Fi(:, Nt_ts + 1));
%             texpE_RHS2 = norm(Fi_testpoint - Fi(:, Nt_ts + 1))/norm(-1i*Hpsi(K, Vthalf, Unew(:, Nt_ts + 1)) +  Fi(:, Nt_ts + 1));
%            texpansion_error3 = norm(Fi_testpoint - Fi(:, Nt_ts + 1))*test_tpoint/norm(Unew(:, Nt_ts + 1));            
%             Fi_testpoint = Cnewton(:, Nt_ts);
%             for ti = (Nt_ts-1):-1:1
%                 Fi_testpoint = Cnewton(:, ti) + (4/Tts)*(test_tpoint - t_ts(ti))*Fi_testpoint;
%             end
%            texpansion_error3 = abs_texpEfactor*norm(Cnewton(:, Nt_ts + 1))*test_tpoint/norm(Unew(:, Nt_ts + 1));
            texpansion_error11 = abs_texpEfactor*norm(Cnewton(:, Nt_ts + 1))*Tts/norm(Unew(:, Nt_ts + 1));
% This may work for even Nt_ts:
            texpansion_error3 = texpEfactor1*norm(Cnewton(:, Nt_ts + 1))/norm(Unew(:, Nt_ts));
% This may work for odd Nt_ts:
            DCnewton = (Cnewton_der(:, 2) - Cnewton_der(:, 1))/(t_integ(Ninteg/2 + 1) - t_integ(Ninteg/2));
            texpansion_error6 = texpEfactor2*norm(DCnewton)/norm(Unew(:, Nt_ts));
            HCnewton = 1i*Hpsi(K, Vts, Cnewton_der(:, 1));
            %[Vkr, Hkr] = createKrop(@(u) Hpsi(K, Vts, u), Lambda(:, Nt_ts + 1), Ncheb);
            texpansion_error8 = texpEfactor2*norm(DCnewton + HCnewton)/norm(Unew(:, Nt_ts));
            texpansion_error9 = texpEfactor2*norm(HCnewton)/norm(Unew(:, Nt_ts));
            % This is a quite good estimation. It doesn't matter much if
            % the Hessenbenrg matrix and lambda come from the previous
            % time-step:
%             DCnewtonKr = Vkr(:, 1:Ncheb)'*DCnewton;
%             Cnewton_derKr = Vkr(:,1:Ncheb)'*Cnewton_der(:, 1);
%             texpansion_error12 = texpEfactor2*norm(DCnewtonKr*norm(DCnewton)/norm(DCnewtonKr) +...
%                 1i*Hkr(1:Ncheb, 1:Ncheb)*Cnewton_derKr*norm(Cnewton_der(:, 1))/norm(Cnewton_derKr))/norm(Unew(:, Nt_ts));
            % This is required, since the projection on the Krylov space
            % causes loss of norm.
            %texpansion_error12 = texpEfactor2*norm(Vkr(:, 1:Ncheb)'*DCnewton + 1i*Hkr(1:Ncheb, 1:Ncheb)*Vkr(:,1:Ncheb)'*Cnewton_der(:, 1))/norm(Unew(:, Nt_ts));
            texpansion_error12 = texpEfactor2*(norm(DCnewton) + Edomain(2)*norm(Cnewton_der(:, 1)))/norm(Unew(:, Nt_ts));
            if tsi>1
                DCnewton1 = (Cnewton(:, Nt_ts + 1) - Cnewton_last)/Tts/(Nt_ts + 1);
                texpansion_error7 = texpEfactor2*norm(DCnewton1)/norm(Unew(:, Nt_ts));
% This seems to be an awful estimation:
%                 Cnewton_lastFi = new_divdif([t_ts, t_ts((Nt_ts - 2):(Nt_ts - 1)) - Tts]*4/Tts,...
%                     Fi_last(:, 1:2) + (Fi(:, 1) - Fi_last(:, 3))*ones(1, 2), diagonal);
%                 DCnewton2 = (Cnewton_lastFi(:, 2) - Cnewton_lastFi(:, 1))/(t_ts(Nt_ts - 1) - t_ts(Nt_ts - 2));
%                 texpansion_error8 = texpEfactor2*norm(DCnewton2)/norm(Unew(:, Nt_ts));
%                texpansion_error8 = texpEfactor1*norm(Cnewton_lastFi(:, 2))/norm(Unew(:, Nt_ts));
            end
            %texpansion_error4 = 4*norm(sum(Cnewton_int*diag(Sn(theta_int, Nt_ts)), 2))*(Tts/Ninteg)/norm(Unew(:, Nt_ts));
            Fi_int_interp = NewtonIpln4b(t(1:Nt_ts), Fi(:, 1:Nt_ts), t_integ);
            texpansion_error5 = norm(sum(Fi_int_interp - Fi_int, 2))*(Tts/Ninteg)/norm(Unew(:, Nt_ts));
            % For a nonlinear problem, Vthalf should be fixed.
            Uinterp = NewtonIpln4b(t(1:Nt_ts), Unew(:, 1:Nt_ts), t(Nt_ts + 1));
            UintE = norm(Uinterp - Unew(:, Nt_ts + 1))/norm(Unew(:, Nt_ts + 1));
%            [tsi, reldif, reldifD, texpansion_error, texpE_RHS, texpansion_error2, texpE_RHS2, UintE]
            HexpE = norm(CchebFts(Ncheb, Nt_ts)*Vcheb(:, Ncheb))/norm(Vcheb*CchebFts(:, Nt_ts));
            Ulast = Unew;
            Ulast_int = Unew_int;
            niter = niter + 1;
        end
        FUerror = Ferror*nFnU;
        for ti = 1:Nt_ts
            Fi_new(:, ti) = -1i*(Vtfun(Unew(:, ti), x, t(ti), varargin{:}) - Vthalf).*Unew(:, ti);
        end
        % If there is an inhomogeneous term, we add it to the fi
        % vectors:
        if there_is_ih
            Fi_new = Fi_new + ihterm(:, 1:Nt_ts);
        end
        conv_error1 = norm(Fi_new(:, Nt_ts) - Fi(:, Nt_ts))*Tts/norm(Unew(:, Nt_ts));
        conv_error2 = norm((Fi_new(:, 1:Nt_ts) - Fi(:, 1:Nt_ts))*wcheb)/norm(Unew(:, Nt_ts));
        %errors(tsi, :) = [reldif, reldifD, texpansion_error, texpE_RHS, texpansion_error2, texpE_RHS2, texpansion_error3, UintE, HexpE, FUerror];
        errors(tsi, :) = [reldif, reldifD, conv_error1, conv_error2, texpansion_error, texpE_RHS, texpansion_error3, texpansion_error5, texpansion_error6,...
            texpansion_error7, texpansion_error8, texpansion_error9, texpansion_error10, texpansion_error11, texpansion_error12, UintE, HexpE, FUerror];
%        texp_dif(:, tsi) = Fi_testpoint - Fi(:, Nt_ts + 1);
        texp_dif(:, tsi) = texpEfactor*Cnewton(:, Nt_ts + 1);
        if niter == Niter
            display('The program has failed to achieve the desired tolerance.')
            % In such a case, change Nts, Nt_ts and/or Ncheb.
        end
        allniter = allniter + niter;
        U(:, tsi+1) = Unew(:, Nt_ts);
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
%        Uguess(:, 2:Nt_ts) = UfromLamb(Lambda, timeMnext, Vcheb, CchebFnext);
        Uguess(:, 2:(Nt_ts + 1)) = UfromLamb(Lambda, timeMnext, Vcheb, CchebFnext);
        Uguess_int(:, 1:Ninteg) = UfromLamb(Lambda, timeMnext_int, Vcheb, CchebFnext_int);
        Cnewton_last = Cnewton(:, Nt_ts + 1);
%         Fi_last = Fi(:, (Nt_ts - 2):Nt_ts);
        if ~isfinite(U(1, tsi + 1))
            % It means the solution diverges.
            % In such a case, change Nts, Nt_ts and/or Ncheb.
            display('Error.');
            return
        end
        if there_is_ih
            ihterm(:, 1) = ihterm(:, Nt_ts);
        end
    end    
    mniter = allniter/Nts;
    matvecs = allniter*(Nt_ts + Ncheb-1);
end

% function [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb)
% % The function computes the coefficients for the Chebyshev expansion of the
% % function of matrix: F(x, t).
% % Note that the coefficients are depndent only on the function, and not on the
% % argument (the operator), or the vector that the function operates on.
% %    Nt_2ts = 2*Nt_ts - 1;
%     Nt_2ts = 2*Nt_ts + 1;
%     CchebF = zeros(Ncheb, Nt_2ts);
%     % For every time point, we consider the function as a different
%     % function of x only, and compute the coefficients.
%     % The coefficients are computed by the function: chebc.
%     % Every column of CchebF represents the coefficients of F(x, t = t_ts(ti)).
%     for ti = 1:Nt_2ts
%         CchebF(:, ti) = chebc(@(x) Ffun(x, t_2ts(ti), Nt_ts, Ncheb), minE, maxE, Ncheb).';      
%     end
% %     CchebFts = CchebF(:, 1:Nt_ts);
% %     CchebFnext = CchebF(:, (Nt_ts + 1):(Nt_2ts));
%     CchebFts = CchebF(:, 1:(Nt_ts + 1));
%     CchebFnext = CchebF(:, (Nt_ts + 2):(Nt_2ts));
% end

function [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb, tol)
    Nt_2ts = length(t_2ts);
    Echeb = cos(((1:Ncheb).'*2 - 1)*pi/(2*Ncheb));
    E = 0.5*(Echeb*(maxE - minE) + maxE + minE);
    F_Et = Ffun(E, t_2ts, Nt_ts, tol);
    CchebF = chebcM(F_Et);
%     CchebFts = CchebF(:, 1:Nt_ts);
%     CchebFnext = CchebF(:, (Nt_ts + 1):(2*Nt_ts - 1));
    CchebFts = CchebF(:, 1:(Nt_ts + 1));
%    CchebFnext = CchebF(:, (Nt_ts + 2):(2*Nt_ts + 1));   
    CchebFnext = CchebF(:, (Nt_ts + 2):Nt_2ts);   
end

function CchebF = FchebC1(t, Nt_ts, minE, maxE, Ncheb, tol)
% The function computes the Chebyshev coefficients of F(x,t), where t serves as a
% parameter. t represents a row vector of time-values.
% The domain of approximation is the eigenvalue domain, [minE, maxE].
    Echeb = cos(((1:Ncheb).'*2 - 1)*pi/(2*Ncheb));
    E = 0.5*(Echeb*(maxE - minE) + maxE + minE);
    F_Et = Ffun(E, t, Nt_ts, tol);
    CchebF = chebcM(F_Et);
end

function timeM = maketM1(t, Nt_ts)
% Computation of the matrix timeM, of time Taylor polynomials.
% timeM(i, j) = t(j)^(i - 1)
    Nt = length(t);
    timeM = zeros(Nt_ts, Nt);
    timeM(1, :) = ones(1, Nt);
    for lambdai = 2:Nt_ts
        timeM(lambdai, :) = t.*timeM(lambdai - 1, :);
    end
end



function result = Ffun(x, t, Nt_ts, tol)
% The function of matrix F(x, t):
    Nt = length(t);
    Nx = length(x);
    minus_ixt = -1i*x*t;
    is_big = factorial(Nt_ts)*eps./abs(minus_ixt.^(Nt_ts)) < tol;
    result = ones(Nx, Nt);
    result(is_big) = exp(minus_ixt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_ixt(is_big);
    end
    is_not_converged = ~is_big;
    term = double(is_not_converged);
    polydeg = 1;
    while max(max(is_not_converged))
        term(is_not_converged) = minus_ixt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
    end
    result = result.*((ones(Nx, 1)*t).^Nt_ts);
end

function [timeMts, timeMnext] = maketM(t_2ts, Nt_ts)
% Computation of the matrix timeM, of time Taylor polynomials.
% timeM(i, j) = t_2ts(j)^(i - 1)
%    Nt_2ts = 2*Nt_ts - 1;
    Nt_2ts = 2*Nt_ts + 1;
    %Nt_2ts = length(t_2ts);
    timeM = zeros(Nt_ts, Nt_2ts);
    timeM(1, :) = ones(1, Nt_2ts);
    for lambdai = 2:Nt_ts
        timeM(lambdai, :) = t_2ts.*timeM(lambdai - 1, :);
    end
%     timeMts = timeM(:, 1:Nt_ts);
%     timeMnext = timeM(:, (Nt_ts + 1):Nt_2ts);    
    timeMts = timeM(:, 1:(Nt_ts + 1));
    timeMnext = timeM(:, (Nt_ts + 2):Nt_2ts);    
end

function Uguess = guess_ts1(ui, Nt_ts)
    Uguess = ui*ones(1, Nt_ts);
end

function [U, nFnU] = UfromLamb(Lambda, timeM, Polyv, CchebF)
    FHlamb = Polyv*CchebF;
%    U = Lambda(:, 1:(end-1))*timeM + Nt_tsPolyv*CchebF;
    U = Lambda(:, 1:(end-1))*timeM + FHlamb;
    if nargout>1
        nFnU = norm(FHlamb(:, end - 1))/norm(U(:, end - 1));
    end
end

function result = Sn(theta, n)
    result = 0.5*sin((n - 1)*theta).*(cos(2*theta) - 1)./sin(theta);
end
