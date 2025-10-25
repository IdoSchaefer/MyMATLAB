function [U mniter] = TDHdiag_tsnf(H0, Vt, ui, x, T, Nts, Nt_ts, tol)
% The program solves Schrodinger equation for a time dependent Hamiltonian.
% H0 is the time-independent Hamiltonian matrix.
% Vt is the time-dependent disturbance. It is a function handle of the form: @(x, t).
% Intended for the case when Vt is a function of the variable of
% ui, and has a known analytical form.
% x is the grid of ui.
% ui is the initial state vector.
% T is the total time. Nts is the number of time steps. Nt_ts is the number
% of Chebychv points within the time step.
% tol is the desired tolerance of the convergence.
% The program is based the "Chebychev propagator with iterative time ordering"
% procedure. The matrices of functions are caculated using a
% diagonalization method.
    Tts = T/Nts;
    dim = length(ui);
    U = zeros(dim, Nts + 1);
    [P, D] = eig(H0);
    U(:, 1) = P\ui;
    eigval = diag(D);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*Tts;
    t_2ts = [t_ts.'; Tts + t_ts(2:Nt_ts).'].';
%     Nt_2ts = 2*Nt_ts - 1;
    expeigval = exp(-1i*eigval*t_2ts);
    [Fts, Fnext] = makeF(eigval, expeigval, dim, Nt_ts, t_2ts);
%    Cc2t = cheb2taylor(Nt_ts-1);
%    Cr2t = r2Taylor(t_ts);
    Cr2t = r2Taylor4(t_ts, Tts);
%    U0 = zeros(dim, Nt_ts+1);
    Unew = zeros(dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    Lambda = zeros(dim, Nt_ts+1);
%    Ctaylor = zeros(dim, Nt_ts);
    Uguess = guess_ts1(ui, expeigval(:, 1:Nt_ts), dim, Nt_ts);
    allniter = 0;
    Niter = 1000;
    for tsi = 1:Nts
        Ulast = Uguess;
        Lambda(:, 1) = U(:, tsi);
%        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Calculation of the inhomogeneous fi vectors:
            for ti = 1:Nt_ts
                Fi(:, ti) = -1i*(P\(Vt(x, (tsi-1)*Tts + t_ts(ti)).*(P*Ulast(:, ti))));
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation in the points t_ts:
            Cnewton = devdif(t_ts*4/Tts, Fi);
%            Ctaylor = Newton2Taylor(t_ts, Cnewton);            
            Ctaylor = zeros(dim, Nt_ts);
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
%                Lambda(:, polyi) = -1i*eigval.*Lambda(:, polyi-1) + Ctaylor(:, polyi-1);
                Lambda(:, polyi) = (-1i*eigval.*Lambda(:, polyi-1) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            % Calculation of the wave function in all the time points
            % within the time step:
            Unew = UfromLamb(Lambda, Fts, dim, t_ts, Nt_ts);
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            fprintf('The program has failed to achieve the desired tolerance.\n')
        end
        allniter = allniter + niter;
        %niter
        U(:, tsi+1) = Unew(:, Nt_ts);
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        Uguess(:, 2:Nt_ts) = UfromLamb(Lambda, Fnext, dim, Tts + t_ts(2:Nt_ts), Nt_ts-1);
        if ~isfinite(U(1, tsi + 1))
            display('Error.');
            return
        end
    end
    U(:, 1) = ui;
    U(:, 2:Nts+1) = P*U(:, 2:Nts+1);
    mniter = allniter/Nts;
end

function [Fts, Fnext] = makeF(eigval, expeigval, dim, Nt_ts, t_2ts) 
    Nt_2ts = 2*Nt_ts - 1;
    F = zeros(dim, Nt_2ts);
    tsmall = t_2ts(min(abs(eigval.'))*t_2ts<Nt_ts+1);
    Ntsmall = length(tsmall);
    tlarge = t_2ts(min(abs(eigval.'))*t_2ts>=Nt_ts+1);
    Ntlarge = Nt_2ts - Ntsmall;
    Mtaylor = -1i*eigval*tsmall;
    for ti = 1:Ntsmall
        Mterm = ones(dim, 1);
%        lastMterm = Mterm;
        F(:, ti) = Mterm;
        polydeg = 1;
        while norm(Mterm)/norm(F(:, ti)) > 1e-16
%        while norm(Mterm)/norm(lastMterm) > 1e-3
%            lastMterm = Mterm;
            Mterm = Mtaylor(:, ti).*Mterm/(polydeg + Nt_ts);
            F(:, ti) = F(:, ti) + Mterm;
            polydeg = polydeg + 1;
        end
    end
    F(:, (Ntsmall + 1):Nt_2ts) = expeigval(:, (Ntsmall + 1):Nt_2ts);
    Mezer = -1i*eigval*tlarge;
    for polyi = 1:Nt_ts
        F(:, (Ntsmall + 1):Nt_2ts) = polyi*(F(:, (Ntsmall + 1):Nt_2ts) - ones(dim, Ntlarge))./Mezer;
    end
    F = F.*(ones(dim, 1)*(t_2ts.^Nt_ts));
    Fts = F(:, 1:Nt_ts);
    Fnext = F(:, (Nt_ts + 1):(Nt_2ts));
end

function Uguess = guess_ts1(ui, expeigval, dim, Nt_ts)
    Uguess = zeros(dim, Nt_ts);
    for ti = 1:Nt_ts
        Uguess(:, ti) = expeigval(:, ti).*ui;
    end
end

function U = UfromLamb(Lambda, F, dim, t, Nt)
    U = zeros(dim, Nt);
    for ti = 1:Nt
        U(:, ti) = F(:, ti).*Lambda(:, Nt+1);
    end
    for polydeg = 0:(Nt-1)
        U = U + Lambda(:, polydeg+1)*t.^polydeg;
    end
end