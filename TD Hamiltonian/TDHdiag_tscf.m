function U = TDHdiag_tscf(H0, Vt, ui, x, T, Nts, Nt_ts, tol)
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
    szx = size(x);
    if szx(1) == 1
        x = x.';
    end
    [P, D] = eig(H0);
    U(:, 1) = P\ui;
    eigval = diag(D);
    tcheb = cos(((1:Nt_ts)*2 - 1)*pi/(2*Nt_ts));
    t_ts = 0.5*(tcheb+1)*Tts;
    allt = [t_ts Tts];
    expeigval = exp(-1i*eigval*allt);
    FNt_ts = expeigval;
    for polydeg = 0:(Nt_ts - 1)
        FNt_ts = FNt_ts - (-1i*eigval*allt).^polydeg/factorial(polydeg);
    end
    vezer = (-1i*eigval).^(-Nt_ts);
    for ti = 1:(Nt_ts + 1)
        FNt_ts(:, ti) = FNt_ts(:, ti).*vezer;
    end
    Cc2t = cheb2taylor(Nt_ts-1);
%    U0 = zeros(dim, Nt_ts+1);
    Ulast = zeros(dim, Nt_ts+1);    
    Unew = zeros(dim, Nt_ts+1);
    Fi = zeros(dim, Nt_ts);
    Lambda = zeros(dim, Nt_ts+1);
    Ctaylor = zeros(dim, Nt_ts);
    Niter = 10000;
    for tsi = 1:Nts
        ui_ts = U(:, tsi);
        for ti = 1:Nt_ts+1
            Ulast(:, ti) = expeigval(:, ti).*ui_ts;
        end
        Lambda(:, 1) = ui_ts;
%        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Calculation of the fi vectors:
            for ti = 1:Nt_ts
                Fi(:, ti) = -1i*(P\(Vt(x, (tsi-1)*Tts + t_ts(ti)).*(P*Ulast(:, ti))));
            end
            % Calculation of the Chebychev coefficients for expansion in
            % time:
            Cchebt = dct(Fi.').'/sqrt(Nt_ts);
            Cchebt(:, 2:Nt_ts) = Cchebt(:, 2:Nt_ts)*sqrt(2);
            %Calculation of the corresponding coefficients for a Taylor
            %expansion:
            Ctaylor(:, Nt_ts) = (2/Tts)^(Nt_ts-1)*Cc2t(Nt_ts, Nt_ts)*Cchebt(:, Nt_ts);
            for polyi = (Nt_ts-1):-1:1
                Ctaylor(:, polyi) = Cc2t(polyi, polyi)*Cchebt(:, polyi);
                for polyj = (polyi+1):Nt_ts
                    Ctaylor(:, polyi) = Ctaylor(:, polyi) + Cc2t(polyj, polyi)*Cchebt(:, polyj) -...
                        (Tts/2)^(polyj-1)/factorial(polyj - polyi)*Ctaylor(:, polyj);
                end
                Ctaylor(:, polyi) = (2/Tts)^(polyi-1)*Ctaylor(:, polyi);
            end
            % Calculation of the lambda vectors:
            for polyi = 2:(Nt_ts+1)
                Lambda(:, polyi) = -1i*eigval.*Lambda(:, polyi-1) + Ctaylor(:, polyi-1);
            end
            % Calculation of the wave function in all the time points
            % within the time step:       
            for ti = 1:Nt_ts+1
                Unew(:, ti) = FNt_ts(:, ti).*Lambda(:, Nt_ts+1);
                for polydeg = 0:(Nt_ts-1)
                    Unew(:, ti) = Unew(:, ti) + allt(ti)^polydeg/factorial(polydeg)*Lambda(:, polydeg+1);
                end
            end
            reldif = norm(Unew(:, Nt_ts+1) - Ulast(:, Nt_ts+1))/norm(Ulast(:, Nt_ts+1));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            fprintf('The program has failed to reach the desired tolerance.\n')
        end
        %niter
%         % Calculation of the wave function in the edge of the time step:
%         U(:, tsi + 1) = FNt_ts(:, ti).*Lambda(:, Nt_ts+1);
%         for polydeg = 0:(Nt_ts-1)
%             U(:, tsi + 1) = U(:, tsi + 1) + t_ts(ti)^polydeg/factorial(polydeg)*Lambda(polydeg+1);
%         end
        U(:, tsi+1) = Unew(:, Nt_ts+1);
    end
    U(:, 1) = ui;
    U(:, 2:Nts+1) = P*U(:, 2:Nts+1);
end