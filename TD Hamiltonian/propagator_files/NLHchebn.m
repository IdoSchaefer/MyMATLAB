function [U mniter matvecs] = NLHchebn(V, Vt, Edomain, ui, xdomain, T, Nts, Nt_ts, Ncheb, tol)
% The program solves Schrodinger equation for a time dependent nonlinear Hamiltonian.
% V is the time-independent linear potential. It's a function handle of the form:
% @(x).
% Vt is the time-dependent nonlinear disturbance. It is a function handle of the form: @(u, x, t).
% Intended for the case when Vt is a function of the variable of
% ui, and has a known analytical form.
% x is the grid of ui.
% ui is the initial state vector.
% T is the total time. Nts is the number of time steps. Nt_ts is the number
% of Chebychv points within the time step.
% tol is the desired tolerance of the convergence.
% The program is based the "Chebychev propagator with iterative time ordering"
% procedure.
    Tts = T/Nts;
    tmidi = round(Nt_ts/2);
    Nx = length(ui);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
    Vvec = V(x);
    minE = Edomain(1);
    maxE = Edomain(2);
    U = zeros(Nx, Nts + 1);
    U(:, 1) = ui;
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*Tts;
    t_2ts = [t_ts.'; Tts + t_ts(2:Nt_ts).'].';
%     Nt_2ts = 2*Nt_ts - 1;
    [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb);
    [timeMts, timeMnext] = maketM(t_2ts, Nt_ts);
%    Cc2t = cheb2taylor(Nt_ts-1);
%    Cr2t = r2Taylor(t_ts);
    Cr2t = r2Taylor4(t_ts, Tts);
%    U0 = zeros(Nx, Nt_ts+1);
    Unew = zeros(Nx, Nt_ts);
    Fi = zeros(Nx, Nt_ts);
    Lambda = zeros(Nx, Nt_ts+1);
%    Ctaylor = zeros(Nx, Nt_ts);
    Uguess = guess_ts1(ui, Nt_ts);
    allniter = 0;
%    matvecs = 0;
    Niter = 1000;
    for tsi = 1:Nts
        Ulast = Uguess;
        Lambda(:, 1) = U(:, tsi);
%        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            Vhalfv = Vt(Ulast(:, tmidi), x, (tsi-1)*Tts + t_ts(tmidi));
            Vts = Vvec + Vhalfv;
            % Calculation of the inhomogeneous fi vectors:
            for ti = 1:Nt_ts
                Fi(:, ti) = -1i*(Vt(Ulast(:, ti), x, (tsi-1)*Tts + t_ts(ti)) - Vhalfv).*Ulast(:, ti);
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation in the points t_ts:
            Cnewton = devdif(t_ts*4/Tts, Fi);
%            Ctaylor = Newton2Taylor(t_ts, Cnewton);            
            Ctaylor = zeros(Nx, Nt_ts);
            Ctaylor(:, 1) = Cnewton(:, 1);
            for Newtoni = 2:Nt_ts
                for Taylori = 1:Newtoni
                    Ctaylor(:, Taylori) = Ctaylor(:, Taylori) + Cnewton(:, Newtoni)*Cr2t(Newtoni, Taylori);
                end
            end
%             % Calculation of the Chebychev coefficients for expansion in
%             % time:
%             Cchebt = dct(Fi.').'/sqrt(Nt_ts);
%             Cchebt(:, 2:Nt_ts) = Cchebt(:, 2:Nt_ts)*sqrt(2);
%             %Calculation of the corresponding coefficients for a Taylor
%             %expansion:
%             Ctaylor(:, Nt_ts) = (2/Tts)^(Nt_ts-1)*Cc2t(Nt_ts, Nt_ts)*Cchebt(:, Nt_ts);
%             for polyi = (Nt_ts-1):-1:1
%                 Ctaylor(:, polyi) = Cc2t(polyi, polyi)*Cchebt(:, polyi);
%                 for polyj = (polyi+1):Nt_ts
%                     Ctaylor(:, polyi) = Ctaylor(:, polyi) + Cc2t(polyj, polyi)*Cchebt(:, polyj) -...
%                         (Tts/2)^(polyj-1)/factorial(polyj - polyi)*Ctaylor(:, polyj);
%                 end
%                 Ctaylor(:, polyi) = (2/Tts)^(polyi-1)*Ctaylor(:, polyi);
%             end
            % Calculation of the Lambda vectors:
            for polyi = 2:(Nt_ts+1)
                Lambda(:, polyi) = (-1i*Hpsi(K, Vts ,Lambda(:, polyi-1)) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            % Calculation of the wave function in all the time points
            % within the time step:
            Vcheb = vchebM(K, Vts, Lambda(:, Nt_ts + 1), Nx, minE, maxE, Ncheb);
            Unew = UfromLamb(Lambda, timeMts, Vcheb, CchebFts);
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            fprintf('The program has failed to achieve the desired tolerance.\n')
        end
        allniter = allniter + niter;
%        matvecs = matvecs + niter*(Nt_ts + Ncheb-1);
        %niter
        U(:, tsi+1) = Unew(:, Nt_ts);
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        Uguess(:, 2:Nt_ts) = UfromLamb(Lambda, timeMnext, Vcheb, CchebFnext);
        if ~isfinite(U(1, tsi + 1))
            display('Error.');
            return
        end
    end    
    mniter = allniter/Nts;
    matvecs = allniter*(Nt_ts + Ncheb-1);
end

function [CchebFts, CchebFnext] = FchebC(t_2ts, Nt_ts, minE, maxE, Ncheb) 
    Nt_2ts = 2*Nt_ts - 1;
    CchebF = zeros(Ncheb, Nt_2ts);
    for ti = 1:Nt_2ts
        CchebF(:, ti) = chebc(@(x) Ffun(x, t_2ts(ti), Nt_ts, Ncheb), minE, maxE, Ncheb).';      
    end
%     tsmall = t_2ts(min(abs(eigval.'))*t_2ts<Nt_ts+1);
%     Ntsmall = length(tsmall);
%     for ti = 1:Ntsmall
%         CchebF(:, ti) = chebc(@(x) Ftaylor(x, t_2ts(ti), Nt_ts), minE, maxE, Ncheb).';      
%     end
%     for ti = (Ntsmall + 1):Nt_2ts
%         CchebF(:, ti) = chebc(@(x) Fdirect(x, t_2ts(ti), Nt_ts), minE, maxE, Ncheb).';
%     end
    CchebFts = CchebF(:, 1:Nt_ts);
    CchebFnext = CchebF(:, (Nt_ts + 1):(Nt_2ts));
end

function result = Ffun(x, t, Nt_ts, Ncheb)
    minus_ixt = -1i*x*t;
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