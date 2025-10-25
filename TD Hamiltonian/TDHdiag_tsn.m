function U = TDHdiag_tsn(H0, Vt, ui, T, Nts, Nt_ts, tol)
% The program solves Schrodinger equation for a time dependent Hamiltonian.
% H0 is the time-independent Hamiltonian matrix.
% Vt is the time-dependent disturbance. It is a 4D tensor; the first 2 indices represent the indices
% of the matrix in a particular time point. The 3'rd index represents the
% time step. The 4'th represents the Chebychev time points within the time
% step, of the Nt_ts order Chebychev polynom.
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
    expeigval = exp(-1i*eigval*t_ts);
% This procedure to obtain F is instable for small t:
%     FNt_ts = expeigval;
%     for polydeg = 0:(Nt_ts-1)
%         FNt_ts = FNt_ts - (-1i*eigval*t_ts).^polydeg/factorial(polydeg);        
%     end
%     FNt_ts = FNt_ts*factorial(Nt_ts);
%     vezer = (-1i*eigval).^(-Nt_ts);
%     for ti = 1:Nt_ts
%         FNt_ts(:, ti) = FNt_ts(:, ti).*vezer;
%     end
% This is another option.
    FNt_ts = zeros(dim, Nt_ts);
    tsmall = t_ts(min(abs(eigval.'))*t_ts<Nt_ts+1);
    Ntsmall = length(tsmall);
    tlarge = t_ts(min(abs(eigval.'))*t_ts>=Nt_ts+1);
    Ntlarge = Nt_ts - Ntsmall;
    Mtaylor = -1i*eigval*tsmall;
    for ti = 1:Ntsmall
        Mterm = ones(dim, 1);
%        lastMterm = Mterm;
        FNt_ts(:, ti) = Mterm;
        polydeg = 1;
        while norm(Mterm)/norm(FNt_ts(:, ti)) > 1e-16
%        while norm(Mterm)/norm(lastMterm) > 1e-3
%            lastMterm = Mterm;
            Mterm = Mtaylor(:, ti).*Mterm/(polydeg + Nt_ts);
            FNt_ts(:, ti) = FNt_ts(:, ti) + Mterm;
            polydeg = polydeg + 1;
        end
    end
    FNt_ts(:, (Ntsmall + 1):Nt_ts) = expeigval(:, (Ntsmall + 1):Nt_ts);
    Mezer = -1i*eigval*tlarge;
    for polyi = 1:Nt_ts
        FNt_ts(:, (Ntsmall + 1):Nt_ts) = polyi*(FNt_ts(:, (Ntsmall + 1):Nt_ts) - ones(dim, Ntlarge))./Mezer;
    end
    FNt_ts = FNt_ts.*(ones(dim, 1)*(t_ts.^Nt_ts));
%%% and another option, instable for small t.
%     Meig = -1i*eigval*ones(1, Nt_ts);
%     Mtaylor = (ones(dim, 1)*t_ts).*Meig;
%     Ftaylor = ones(dim, Nt_ts);
%     Mterm = ones(dim, Nt_ts);
%     for polydeg = 1:Nt_ts-1
%         Mterm = Mtaylor.*Mterm/polydeg;
%         Ftaylor = Ftaylor + Mterm;
%     end
%     FNt_ts = FNt_ts - Ftaylor;
%     for polyi = 1:Nt_ts
%         FNt_ts = FNt_ts*polyi./Meig;
% %        Ftaylor = Ftaylor*polyi./Meig;
%     end
%    Cc2t = cheb2taylor(Nt_ts-1);
%    Cr2t = r2Taylor(t_ts);
    Cr2t = r2Taylor1(t_ts);
%    U0 = zeros(dim, Nt_ts+1);
    Ulast = zeros(dim, Nt_ts);
    Unew = zeros(dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    Lambda = zeros(dim, Nt_ts+1);
%    Ctaylor = zeros(dim, Nt_ts);
    Niter = 1000;
    for tsi = 1:Nts
        ui_ts = U(:, tsi);
        for ti = 1:Nt_ts
            Ulast(:, ti) = expeigval(:, ti).*ui_ts;
        end
        Lambda(:, 1) = ui_ts;
%        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Calculation of the inhomogeneous fi vectors:
            for ti = 1:Nt_ts
                Fi(:, ti) = -1i*(P\Vt(:, :, tsi, ti)*P*Ulast(:, ti));
            end
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation in the points t_ts:
            Cnewton = devdif(t_ts, Fi);
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
            % Calculation of the lambda vectors:
            for polyi = 2:(Nt_ts+1)
%                Lambda(:, polyi) = -1i*eigval.*Lambda(:, polyi-1) + Ctaylor(:, polyi-1);
                Lambda(:, polyi) = (-1i*eigval.*Lambda(:, polyi-1) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            % Calculation of the wave function in all the time points
            % within the time step:
            for ti = 1:Nt_ts
                Unew(:, ti) = FNt_ts(:, ti).*Lambda(:, Nt_ts+1);
            end
            for polydeg = 0:(Nt_ts-1)
                Unew = Unew + Lambda(:, polydeg+1)*t_ts.^polydeg;
            end
%            for ti = 1:Nt_ts
%                 Unew(:, ti) = FNt_ts(:, ti).*Lambda(:, Nt_ts+1);
%                 for polydeg = 0:(Nt_ts-1)
%                     Unew(:, ti) = Unew(:, ti) + t_ts(ti)^polydeg/factorial(polydeg)*Lambda(:, polydeg+1);
%                 end
%             end
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            fprintf('The program has failed to reach the desired tolerance.\n')
        end
        %niter
        U(:, tsi+1) = Unew(:, Nt_ts);
    end
    U(:, 1) = ui;
    U(:, 2:Nts+1) = P*U(:, 2:Nts+1);
end