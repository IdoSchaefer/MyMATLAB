function [U, mniter, matvecs, history] = TDHxp1(K, V, Vtfun, ihfun, Edomain, ui, x, tgrid, Nts, Nt_ts, Ncheb, tol, display_mode, varargin)
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
% tgrid: the time grid of the desired output solution, U. Should be ordered in an
% increasing order for a forward propagation, and a decreasing order for a
% backward propagation.
% Nts: the number of time steps of the propagation. 
% Nt_ts: the number of interior Chebyshev time points in the time step, used during
% the computational process.
% Ncheb: the number of terms used to expand the function of matrix in the
% Chebyshev expansion.
% tol: the desired tolerance of the convergence.
% display_mode: A boolean variable. true means that warnings are displayed
% during the propagation. false means that warnings are displayed only
% before and after the propagation.
% Output:
% U: contains the solution at the time points specified in tgrid, in separate columns.
% mniter: the mean number of iteration for a time step; should be close to
% 1 for ideal efficiency of the algorithm.
% matvecs: the number of Hamiltonian operations.
% history: Structure which contains the history of propagation:
%   history.t: The time grid of propagation
%   history.U: The solution at history.t
%   history.Ferror: The estimation for the error of the computation of the
%   function of matrix
%   history.texp_error:  A row vector which contains the estimation for the error of the time expansion
%   in each time-step.
    if nargin<13
        display_mode = true;
    end
    % In order to detect if the propagation is a forward or a backward
    % propagation:
    direction = sign(tgrid(2) - tgrid(1));
    Nt = length(tgrid);
    tinit = tgrid(1);
    tf = tgrid(Nt);
    % The length of the time interval of the whole propagation(can be negative):
    T = tf - tinit;
    % The length of the time step interval:
    Tts = T/Nts;
    % The index of the middle term in the time step:
    tmidi = round(Nt_ts/2);
    minE = Edomain(1);
    maxE = Edomain(2);
    Nx = length(ui);
    U = zeros(Nx, Nt);
    U(:, 1) = ui;
    % The Chebyshev points for expansion in time, in the domain in which the Chebyshev expansion
    % is defined: [-1 1]
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The Chebyshev points for expansion in time, in the domain of the time
    % variable (changing of variables):
    t_ts = 0.5*(tcheb+1)*Tts;
    test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation:
    t_2ts = [t_ts, test_tpoint, Tts + t_ts(2:Nt_ts), Tts + test_tpoint];
    if nargout>3
        history.U = zeros(Nx, Nts*(Nt_ts - 1) + 1);
        history.U(:, 1) = ui;
        history.t = [kron(0:Tts:(T - Tts), ones(1, Nt_ts - 1)), T] + [kron(ones(1, Nts), t_ts(1:(Nt_ts - 1))), 0];
        history.texp_error = zeros(1, Nts);
    end
    % Computing the coefficients for the Chebyshev expansion of the
    % function of matrix, in all the interior time points.
    % CchebFts contains the coefficients of the current time step, and
    % CchebFnext contains the coefficients of the next one.
    CchebFcomp = FchebC(t_2ts(2:(2*Nt_ts + 1)), Nt_ts, minE, maxE, Ncheb, tol);
    CchebFts = CchebFcomp(:, 1:Nt_ts);
    CchebFnext = CchebFcomp(:, (Nt_ts + 1):(2*Nt_ts));   
    % Error estimation for the computation of the function of matrix:
    EdomainL = maxE - minE;
    xtest = ((minE + EdomainL/(Nx + 1)):EdomainL/(Nx + 1):(maxE - EdomainL/(Nx + 1))).';
    Fxtest = Ffun(xtest, Tts, Nt_ts, tol);
    Ferror = max(abs(chebc2result(CchebFts(:, Nt_ts - 1), Edomain, xtest) - Fxtest)./abs(Fxtest));
    if Ferror>tol
        fprintf('Warning: The estimated error of the computation of the function of matrix (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', Ferror)
    end
    if Ferror>1e-5
        fprintf('Warning: The estimated error of the computation of the function of matrix (%d) is larger than 1e-5.\nInstability in the propagation process may occur.\n', Ferror)
    end
    if nargout>3
        history.Ferror = Ferror;
    end
    % Computing the matrix of the Taylor polynomials in time.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    timeMcomp = maketM(t_2ts(2:(2*Nt_ts + 1)), Nt_ts);
    timeMts = timeMcomp(:, 1:Nt_ts);
    timeMnext = timeMcomp(:, (Nt_ts + 1):(2*Nt_ts));    
    % Computing the coefficients of the transformation from the Newton
    % interpolation polynomial terms, to a Taylor like form:
    Cr2t = r2Taylor4(t_ts, Tts);
    % The solution in the interior points in the time step.
    % Every column represents an interior time point in the time step:
    %Unew = zeros(Nx, Nt_ts + 1);
    % The "inhomogeneos" vectors, including the term that is realy u
    % dependent (appears in all the cases, not only for the inhomogeneous Schrodinger equation):
    Fi = zeros(Nx, Nt_ts + 1);
    % The Lambda vectors are defined in a recursive way, and contain information about
    % the time dependence of the Fi vectors:
    Lambda = zeros(Nx, Nt_ts + 1);
    % If there is no inhomogeneous term in the equation, ihfun is empty, and there_is_ih = false.
    % If there is, there_is_ih = true.
    there_is_ih = ~isempty(ihfun);
    if there_is_ih
        ihterm = zeros(Nx, Nt_ts);
        ihterm(:, 1) = ihfun(0, varargin{:});
    end
    % The 0'th order approximation is the first guess, for the first time step.
    % Every column represents an interior time point in the time step:
    Uguess = guess_ts1(ui, Nt_ts + 1);
    allniter = 0;
    Niter = 20;
    tgrid_lowi = 2;
    tgrid_upi = 1;
    max_texp_error = 0;
    for tsi = 1:Nts
        % The time, represented by the interior time points:
        t = tinit + Tts*(tsi - 1) + [t_ts, test_tpoint];
        % The first guess for the iterative process, for the convergence of the u
        % values. Every column represents an interior time point in the time step:
        Ulast = Uguess;
        Unew(:, 1) = Ulast(:, 1);
        Lambda(:, 1) = Ulast(:, 1);
        if there_is_ih
            % Computing the inhomogeneous term:
            for ti = 2:(Nt_ts + 1)
                ihterm(:, ti) = ihfun(t(ti), varargin{:});
            end
        end
        % Starting an iterative process, until convergence:
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Vthalf is the x vector representing Vt(t = t(tmidi)): 
            Vthalf = Vtfun(Ulast(:, tmidi), x, t(tmidi), varargin{:});
            % Vts is the x vector representing the total potential energy 
            % operator at t = t(tmidi):
            Vts = V + Vthalf;
            % Calculation of the inhomogeneous fi vectors:
            for ti = 1:(Nt_ts + 1)
                Fi(:, ti) = -1i*(Vtfun(Ulast(:, ti), x, t(ti), varargin{:}) - Vthalf).*Ulast(:, ti);
            end
            % If there is an inhomogeneous term, we add it to the fi
            % vectors:
            if there_is_ih
                Fi = Fi + ihterm;
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation in the points t_ts.
            % The divided differences are computed by the function divdif.
            % For numerical stability, we have to transform the time points
            % in the time step, to points in an interval of length 4:
            Cnewton = divdif([t_ts, test_tpoint]*4/Tts, Fi);
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
            Unew(:, 2:(Nt_ts + 1)) = UfromLamb(Lambda, timeMts, Vcheb, CchebFts);
            % Error estimation for the time expansion.
            % A Newton interpolation of Fi in the test time-point:
            Fi_testpoint = Cnewton(:, Nt_ts);
            for ti = (Nt_ts-1):-1:1
                Fi_testpoint = Cnewton(:, ti) + (4/Tts)*(test_tpoint - t_ts(ti))*Fi_testpoint;
            end
            % The error estimation of the time-expansion:
            texpansion_error = norm(Fi_testpoint - Fi(:, Nt_ts + 1))*Tts/norm(Unew(:, Nt_ts + 1));
            if display_mode && texpansion_error>tol 
                fprintf('Warning: The estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate (in time step No. %d).\n', texpansion_error, tsi)
            end
            % To check the convergence of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            display('Warning: The program has failed to achieve the desired tolerance in the iterative process.')
            % In such a case, change Nts, Nt_ts and/or Ncheb.
        end
        if max_texp_error<texpansion_error
            max_texp_error = texpansion_error;
        end
        if nargout>3
            history.texp_error(tsi) = texpansion_error;
        end
        allniter = allniter + niter;
        % Computation of the solution at the tgrid points.
        % Finding the indices of the points within the time step:
        while tgrid_upi<Nt && (t(Nt_ts) - tgrid(tgrid_upi + 1))*direction>=abs(t(Nt_ts))*eps*10
            tgrid_upi = tgrid_upi + 1;
        end
        if tgrid_lowi<=tgrid_upi
            CchebFout = FchebC(tgrid(tgrid_lowi:tgrid_upi) - t(1), Nt_ts, minE, maxE, Ncheb, tol);
            timeMout = maketM(tgrid(tgrid_lowi:tgrid_upi) - t(1), Nt_ts);
            U(:, tgrid_lowi:tgrid_upi) = UfromLamb(Lambda, timeMout, Vcheb, CchebFout);
            tgrid_lowi = tgrid_upi + 1;
        end
        if abs(t(Nt_ts) - tgrid(tgrid_upi + 1))<abs(t(Nt_ts))*eps*10
            tgrid_upi = tgrid_upi + 1;
            U(:, tgrid_upi) = Unew(:, Nt_ts);
            tgrid_lowi = tgrid_upi + 1;
        end
        if nargout>3
            history.U(:, ((tsi - 1)*(Nt_ts - 1) + 2):(tsi*(Nt_ts - 1) + 1)) = Unew(:, 2:Nt_ts);
        end
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        Uguess(:, 2:(Nt_ts + 1)) = UfromLamb(Lambda, timeMnext, Vcheb, CchebFnext);
        if ~isfinite(Unew(:, Nt_ts))
            % It means the solution diverges.
            % In such a case, change Nts, Nt_ts and/or Ncheb.
            display('Error: The solution diverges.');
            return
        end
        if there_is_ih
            ihterm(:, 1) = ihterm(:, Nt_ts);
        end
    end
    if max_texp_error>tol
        fprintf('\nWarning: The maximal estimated error of the time expansion (%d) is larger than the requested tolerance.\nThe solution might be inaccurate.\n', max_texp_error)
    end
    mniter = allniter/Nts;
    matvecs = allniter*(Nt_ts + Ncheb-1);
end

function CchebF = FchebC(t, Nt_ts, minE, maxE, Ncheb, tol)
% The function computes the Chebyshev coefficients of F(x,t), where t serves as a
% parameter. t represents a row vector of time-values.
% The domain of approximation is the eigenvalue domain, [minE, maxE].
    Echeb = cos(((1:Ncheb).'*2 - 1)*pi/(2*Ncheb));
    E = 0.5*(Echeb*(maxE - minE) + maxE + minE);
    F_Et = Ffun(E, t, Nt_ts, tol);
    CchebF = chebcM(F_Et);
end

function result = Ffun(x, t, Nt_ts, tol)
% The function of matrix F(x, t):
    Nt = length(t);
    Nx = length(x);
    minus_ixt = -1i*x*t;
    % Condition for estimating if F(x, t) should be computed directly or by
    % a "tail" of a Taylor expansion:
    is_big = factorial(Nt_ts)*eps./abs(minus_ixt.^(Nt_ts)) < tol;
    result = ones(Nx, Nt);
    % A direct computation for large arguments:
    result(is_big) = exp(minus_ixt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_ixt(is_big);
    end
    % Computation by a Taylor form for small arguments:
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

function timeM = maketM(t, Nt_ts)
% Computation of the matrix timeM, of time Taylor polynomials.
% timeM(i, j) = t(j)^(i - 1)
    Nt = length(t);
    timeM = zeros(Nt_ts, Nt);
    timeM(1, :) = ones(1, Nt);
    for lambdai = 2:Nt_ts
        timeM(lambdai, :) = t.*timeM(lambdai - 1, :);
    end
end

function Uguess = guess_ts1(ui, Nt_ts)
    Uguess = ui*ones(1, Nt_ts);
end

function U = UfromLamb(Lambda, timeM, Polyv, CchebF)
    U = Lambda(:, 1:(end-1))*timeM + Polyv*CchebF; 
end