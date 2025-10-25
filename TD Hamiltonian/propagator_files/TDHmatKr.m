function [U mniter matvecs] = TDHmatKr(H0, Vtfun, ihfun, ui, tdomain, Nts, Nt_ts, Nkr, tol, varargin)
% The program solves Schrodinger equation for a time dependent Hamiltonian.
% It uses the Hillel Tal-Ezer's propagator.
% The Hamiltonian can be nonlinear (u dependent).
% There is an option to include an inhomogeneos term.
% Edomain: the domain of the eigenenergies of the hamiltonian: [minE maxE]
% ui: the initial state vector.
% tdomain: the time interval of the propagation: [ti tf]
% Nts: the number of time steps. 
% Nt_ts: the number of interior Chebyshev time points in the time step, used during
% the computational process.
% Nkr: the number of vectors in the Krylov space, for coputation of the
% function of operator.
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
    Nx = length(ui);
    U = zeros(Nx, Nts + 1);
    U(:, 1) = ui;
    v0kr = zeros(Nkr, 1);
    v0Dbase = zeros(Nkr, 1);
    % The Chebyshev points for expansion in time, in the domain in which the Chebyshev expansion
    % is defined: [-1 1]
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The Chebyshev points for expansion in time, in the domain of the time
    % variable (changing of variables):
    t_ts = 0.5*(tcheb+1)*Tts;
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation:
    t_2ts = [t_ts.'; Tts + t_ts(2:Nt_ts).'].';
    Vkr = zeros(Nx, Nkr + 1);
    Hkr = zeros(Nkr + 1, Nkr);
    P = zeros(Nkr, Nkr);
    D = zeros(Nkr, Nkr);
    % Computing the matrix of the Taylor polynomials in time.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    [timeMts, timeMnext] = maketM(t_2ts(2:end), Nt_ts);
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
        Unew(:, 1) = Ulast(:, 1);
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
            % Vthalf is the matrix representing Vt(t = t(tmidi)): 
            Vthalf = Vtfun(Ulast(:, tmidi), t(tmidi), varargin{:});
            % Hts is the matrix representing Ht(t = t(tmidi)):
            Hts = H0 + Vthalf;
            % Calculation of the inhomogeneous fi vectors:
            for ti = 1:Nt_ts
                Fi(:, ti) = -1i*(Vtfun(Ulast(:, ti), t(ti), varargin{:}) - Vthalf)*Ulast(:, ti);
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
                Lambda(:, polyi) = (-1i*Hts*Lambda(:, polyi-1) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            eigE = get_eigE;
            % Calculation of the wave function in all the time points
            % within the time step:
            Unew(:, 2:Nt_ts) = UfromLamb(timeMts, t_ts(2:Nt_ts));
            % To check the convergence of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            display('The program has failed to achieve the desired tolerance.')
            % In such a case, change Nts, Nt_ts and/or Nkr.
        end
        allniter = allniter + niter;
        U(:, tsi+1) = Unew(:, Nt_ts);
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        Uguess(:, 2:Nt_ts) = UfromLamb(timeMnext, t_2ts((Nt_ts + 1):(2*Nt_ts - 1)));
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
    matvecs = allniter*(Nt_ts + Nkr);
    
    %%% Nested functions: %%%
    
    function eigE = get_eigE
        createHkr(Lambda(:, Nt_ts + 1));
        [P, D] = eig(Hkr(1:Nkr, :));
        v0kr(1) = norm(Lambda(:, Nt_ts + 1));
        v0Dbase = P\v0kr;
        eigE = diag(D);
    end

    function U = UfromLamb(timeM, t_vec)
        U = Lambda(:, 1:(end-1))*timeM + Vkr(:, 1:Nkr)*P*spdiags(v0Dbase, 0, Nkr, Nkr)*Ffun(eigE, t_vec, Nt_ts, Nkr); 
    end

    function createHkr(v0)
    % The function creates the orthonormalized Krylov space and the
    % Hessenberg matrix Hkr.
        Vkr(:, 1) = v0/norm(v0);
        for vj = 1:Nkr
            Vkr(:, vj+1) = Hts*Vkr(:, vj);
            for vi = 1:vj
                Hkr(vi, vj) = Vkr(:, vi)'*Vkr(:, vj+1);
                Vkr(:, vj+1) = Vkr(:, vj+1) - Hkr(vi, vj)*Vkr(:, vi);
            end
            Hkr(vj+1, vj) = norm(Vkr(:, vj+1));
            Vkr(:, vj+1) = Vkr(:, vj+1)/Hkr(vj+1, vj);
        end
    end

end

%%% Sub functions: %%%

function result = Ffun(x, t, Nt_ts, Nkr)
% The function of matrix F(x, t):
    minus_ixt = -1i*x*t;
    % If the argument of the exponent, -ixt, is small, the function
    % shouldn't be computed directly, because of numerical instability. We
    % use the "tail" of a Taylor expansion instead:
    is_small = abs(minus_ixt) < Nt_ts + 1;
    is_big = ~is_small;
    term = double(is_small);
    result = double(is_small);
    polydeg = 1;
    while max(max(is_small))
        term(is_small) = minus_ixt(is_small).*term(is_small)/(polydeg + Nt_ts);
        result(is_small) = result(is_small) + term(is_small);
        polydeg = polydeg + 1;
        is_small(is_small) = abs(term(is_small)./result(is_small)) > eps;
    end
    % If -ixt is large enough, the function is computed directly:
    result(is_big) = exp(minus_ixt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_ixt(is_big);
    end
    result = result.*((ones(Nkr, 1)*t).^Nt_ts);
end

function [timeMts, timeMnext] = maketM(t_2ts, Nt_ts)
% Computation of the matrix timeM, of time Taylor polynomials.
% timeM(i, j) = t_2ts(j)^(i - 1)
    Nt_2ts = 2*Nt_ts - 2;
    timeM = zeros(Nt_ts, Nt_2ts);
    timeM(1, :) = ones(1, Nt_2ts);
    for lambdai = 2:Nt_ts
        timeM(lambdai, :) = t_2ts.*timeM(lambdai - 1, :);
    end
    timeMts = timeM(:, 1:(Nt_ts - 1));
    timeMnext = timeM(:, Nt_ts:Nt_2ts);    
end

function Uguess = guess_ts1(ui, Nt_ts)
    Uguess = ui*ones(1, Nt_ts);
end