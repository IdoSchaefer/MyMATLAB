function U = TDHcheb_tscf1(H0, minE0, maxE0, Vt, ui, x, T, Nts, Nt_ts, Ncheb, tol)
% The program solves Schrodinger equation for a time dependent Hamiltonian.
% H0 is the time-independent Hamiltonian matrix. minE0 and maxE0 are the
% minimal and maximal eigenenergies of H0.
% Vt is the time-dependent disturbance. It is a function handle of the form: @(x, t).
% Intended for the case when Vt is a function of the variable of
% ui, and has a known analytical form.
% x is the grid of ui.
% ui is the initial state vector.
% T is the total time. Nts is the number of time steps. Nt_ts is the number
% of Chebychv points within the time step.
% Ncheb is the number of Chebychev polynomials to expand the function of
% H0.
% tol is the desired tolerance of the convergence.
% The program is based the "Chebychev propagator with iterative time ordering"
% procedure.
    Tts = T/Nts;
    tmidi = round(Nt_ts/2);
    dim = length(ui);
    I = eye(dim);
    szx = size(x);
    if szx(1) == 1
        x = x.';
    end
    U = zeros(dim, Nts + 1);
    U(:, 1) = ui;
    tcheb = cos(((1:Nt_ts)*2 - 1)*pi/(2*Nt_ts));
    t_ts = 0.5*(tcheb+1)*Tts;
    allt = [t_ts Tts];
    % Computation of the Chebychev coefficients for the expansion of the
    % functions of Hts:
    Cchebexp = zeros(Ncheb, Nt_ts+1);
    CchebF = zeros(Ncheb, Nt_ts+1);
    for ti = 1:(Nt_ts+1)
        fexp = @(x) exp(-1i*x*allt(ti));
        Cchebexp(:, ti) = chebc(fexp, minE0, maxE0, Ncheb).';
        fF = fexp;
        for polydeg = 0:(Nt_ts-1)
            fF = @(x) fF(x) - (-1i*x*allt(ti)).^polydeg/factorial(polydeg);
        end
        fF = @(x) (-1i*x).^(-Nt_ts).*fF(x);
        CchebF(:, ti) = chebc(fF, minE0, maxE0, Ncheb).';        
    end
%    H0cheb = (2*H0 - (minE0 + maxE0)*I)/(maxE0 - minE0);    
    Cc2t = cheb2taylor(Nt_ts-1);
%    U0 = zeros(dim, Nt_ts+1);
%    Ulast = zeros(dim, Nt_ts+1);    
    Unew = zeros(dim, Nt_ts+1);
    Fi = zeros(dim, Nt_ts);
    Lambda = zeros(dim, Nt_ts+1);
    Ctaylor = zeros(dim, Nt_ts);
    Niter = 10000;
    for tsi = 1:Nts
        ui_ts = U(:, tsi);
        Vtts = @(x, t) Vt(x, t) - Vt(x, (tsi-1)*Tts + t_ts(tmidi));
        Hts = H0 + diag(Vt(x, (tsi-1)*Tts + t_ts(tmidi)));
        Htscheb = (2*Hts - (minE0 + maxE0)*I)/(maxE0 - minE0);    
% Computation of the first guess: the solution within the time step, for H = Hts:
        v1 = ui_ts;
        v2 = Htscheb*v1;
        Ulast = v1*Cchebexp(1,:) + v2*Cchebexp(2, :);
        for k = 3:Ncheb
            vk = 2*Htscheb*v2 - v1;
            Ulast = Ulast + vk*Cchebexp(k, :);
            v1 = v2;
            v2 = vk;
        end
        Lambda(:, 1) = ui_ts;
%        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Calculation of the inhomogeneous fi vectors:
            for ti = 1:Nt_ts
                Fi(:, ti) = -1i*Vtts(x, (tsi-1)*Tts + t_ts(ti)).*Ulast(:, ti);
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
                Lambda(:, polyi) = -1i*Hts*Lambda(:, polyi-1) + Ctaylor(:, polyi-1);
            end
            % Calculation of the wave function in all the time points
            % within the time step:
            % Caculation of F(t)*Lambda(:, Nt_ts+1) in all the time points:
            v1 = Lambda(:, Nt_ts+1);
            v2 = Htscheb*v1;
            Unew = v1*CchebF(1,:) + v2*CchebF(2, :);
            for k = 3:Ncheb
                vk = 2*Htscheb*v2 - v1;
                Unew = Unew + vk*CchebF(k, :);
                v1 = v2;
                v2 = vk;
            end
            for ti = 1:(Nt_ts+1)
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
        U(:, tsi+1) = Unew(:, Nt_ts+1);
    end
    niter
end