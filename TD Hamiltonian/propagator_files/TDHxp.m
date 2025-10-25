function [U mniter matvecs] = TDHxp(K, V, Vtfun, ihfun, Edomain, ui, x, tdomain, Nts, Nt_ts, Ncheb, tol, varargin)
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
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation:
    t_2ts = [t_ts.'; Tts + t_ts(2:Nt_ts).'].';
    % Computing the coefficients for the Chebyshev expansion of the
    % function of matrix, in all the interior time points.
    % CchebFts contains the coefficients of the current time step, and
    % CchebFnext contains the coefficients of the next one.
    [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb);
    % Computing the matrix of the Taylor polynomials in time.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    [timeMts, timeMnext] = maketM(t_2ts, Nt_ts);
    % Computing the coefficients of the transformation from the Newton
    % interpolation polynomial terms, to a Taylor like form:
    Cr2t = r2Taylor4(t_ts, Tts);
    % The solution in the interior points in the time step.
    % Every column represents an interior time point in the time step:
    Unew = zeros(Nx, Nt_ts);
    % The "inhomogeneos" vectors, including the term that is realy u
    % dependent (appears in all the cases, not only for the inhomogeneous Schrodinger equation):
    Fi = zeros(Nx, Nt_ts);
    % The Lambda vectors are defined in a recursive way, and contain information about
    % the time dependence of the Fi vectors:
    Lambda = zeros(Nx, Nt_ts+1);
    % If there is no inhomogeneous term in the equation, ihfun = [], and the length is 0 (false).
    % If there is, the length is 1 (true).
    there_is_ih = length(ihfun);
    if there_is_ih
        ihterm = zeros(Nx, Nt_ts);
        ihterm(:, 1) = ihfun(0, varargin{:});
    end
    % The 0'th order approximation is the first guess, for the first time step.
    % Every column represents an interior time point in the time step:
    Uguess = guess_ts1(ui, Nt_ts);
    allniter = 0;
    Niter = 100;
    for tsi = 1:Nts
        % The time, represented by the interior time points:
        t = tinit + Tts*(tsi - 1) + t_ts;
        % The first guess for the iterative process, for the convergence of the u
        % values. Every column represents an interior time point in the time step:
        Ulast = Uguess;
        Lambda(:, 1) = U(:, tsi);
        if there_is_ih
            % Computing the inhomogeneous term:
            for ti = 2:Nt_ts
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
            for ti = 1:Nt_ts
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
            % The divided differences are computed by the function devdif.
            % For numerical stability, we have to transform the time points
            % in the time step, to points in an interval of length 4:
            Cnewton = devdif(t_ts*4/Tts, Fi);
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
            Unew = UfromLamb(Lambda, timeMts, Vcheb, CchebFts);
            % To check the convergence of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            display('The program has failed to achieve the desired tolerance.')
            % In such a case, change Nts, Nt_ts and/or Ncheb.
        end
        allniter = allniter + niter;
        U(:, tsi+1) = Unew(:, Nt_ts);
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        Uguess(:, 2:Nt_ts) = UfromLamb(Lambda, timeMnext, Vcheb, CchebFnext);
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

function [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb)
% The function computes the coefficients for the Chebyshev expansion of the
% function of matrix: F(x, t).
% Note that the coefficients are depndent only on the function, and not on the
% argument (the operator), or the vector that the function operates on.
    Nt_2ts = 2*Nt_ts - 1;
    CchebF = zeros(Ncheb, Nt_2ts);
    % For every time point, we consider the function as a different
    % function of x only, and compute the coefficients.
    % The coefficients are computed by the function: chebc.
    % Every column of CchebF represents the coefficients of F(x, t = t_ts(ti)).
    for ti = 1:Nt_2ts
        CchebF(:, ti) = chebc(@(x) Ffun(x, t_2ts(ti), Nt_ts, Ncheb), minE, maxE, Ncheb).';      
    end
    CchebFts = CchebF(:, 1:Nt_ts);
    CchebFnext = CchebF(:, (Nt_ts + 1):(Nt_2ts));
end

function result = Ffun(x, t, Nt_ts, Ncheb)
% The function of matrix F(x, t):
    minus_ixt = -1i*x*t;
    % If the argument of the exponent, -ixt, is small, the function
    % shouldn't be computed directly, because of numerical instability. We
    % use the "tail" of a Taylor expansion instead:
    if min(abs(minus_ixt)) < (Nt_ts + 1)
        term = ones(1, Ncheb);
        result = ones(1, Ncheb);
        polydeg = 1;
        minus_ixt = -1i*x*t;
        while abs(term/result) > eps
            term = minus_ixt.*term/(polydeg + Nt_ts);
            result = result + term;
            polydeg = polydeg + 1;
        end
    % If -ixt is large enough, the function is computed directly:
    else
        minus_ixt = -1i*x*t;
        result = exp(minus_ixt);
        for polyi = 1:Nt_ts
            result = polyi*(result - 1)./minus_ixt;
        end
    end
    result = result*t^Nt_ts;
end

function [timeMts, timeMnext] = maketM(t_2ts, Nt_ts)
% Computation of the matrix timeM, of time Taylor polynomials.
% timeM(i, j) = t_2ts(j)^(i - 1)
    Nt_2ts = 2*Nt_ts - 1;
    timeM = zeros(Nt_ts, Nt_2ts);
    timeM(1, :) = ones(1, Nt_2ts);
    for lambdai = 2:Nt_ts
        timeM(lambdai, :) = t_2ts.*timeM(lambdai - 1, :);
    end
    timeMts = timeM(:, 1:Nt_ts);
    timeMnext = timeM(:, (Nt_ts + 1):Nt_2ts);    
end

function Uguess = guess_ts1(ui, Nt_ts)
    Uguess = ui*ones(1, Nt_ts);
end

function U = UfromLamb(Lambda, timeM, Polyv, CchebF)
    U = Lambda(:, 1:(end-1))*timeM + Polyv*CchebF; 
end