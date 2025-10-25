function U = TDHdiag_ts2(H0, Vt, ui, T, Nt, Nt_ts, tol)
% Vt is a matrix with Nt columns.Intended for the case when Vt is a function of the variable of
% ui; every column represents a time point.
    dt = T/Nt;
    Tts = Nt_ts*dt;
    Nts = Nt/Nt_ts;
    dim = length(ui);
    U = zeros(dim, Nt + 1);
    [P, D] = eig(H0);
    U(:, 1) = P\ui;
    eigval = diag(D);
    expeigval = exp(-1i*eigval*(dt:dt:Tts));
    U0 = zeros(dim, Nt_ts);
    Niter = 10000;
    for tsi = 1:Nts        
        ui_ts = U(:, (tsi-1)*Nt_ts + 1);
        for ti = 1:Nt_ts
            U0(:, ti) = expeigval(:, ti).*ui_ts;
        end
        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            Unew = U0;
            for taui = 1:Nt_ts
                vtau = -1i*(P\(Vt(:,  (tsi - 1)*Nt_ts + taui).*(P*Ulast(:, taui))))*dt;
                % Adding the term of ti = taui to the sum, which is simply vtau:
                Unew(:, taui) = Unew(:, taui) + vtau;
                for ti = (taui+1):Nt_ts
                    Unew(:, ti) = Unew(:, ti) + expeigval(:, ti - taui).*vtau;
                end
            end
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            fprintf('The program has failed to reach the desired tolerance.\n')
        end
        %niter
        U(:, (tsi-1)*Nt_ts + 1 + (1:Nt_ts)) = Unew;
    end
    U(:, 1) = ui;
    U(:, 2:Nt+1) = P*U(:, 2:Nt+1);
end