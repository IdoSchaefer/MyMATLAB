function U = TDHdiag(H0, Vt, ui, T, Nt, tol)
    dt = T/Nt;
    dim = length(ui);
    U = zeros(dim, Nt + 1);
    U(:, 1) = ui;
    U0 = zeros(dim, Nt);
    [P, D] = eig(H0);
    uiDbase = P\ui;
    eigval = diag(D);
    expeigval = exp(-1i*eigval*(dt:dt:T));
    for ti = 1:Nt
        U0(:, ti) = expeigval(:, ti).*uiDbase;
    end
%    expeigval = conj(expeigval);
    Ulast = U0; 
    reldif = tol + 1;
    niter = 0;
    Niter = 10000;
    while (reldif>tol && niter<Niter)
        Unew = U0;
        for taui = 1:Nt
            vtau = -1i*(P\Vt(:, :, taui)*P*Ulast(:, taui))*dt;
            % Adding the term of ti = taui to the sum, which is simply vtau:
            Unew(:, taui) = Unew(:, taui) + vtau;
            for ti = (taui+1):Nt
                Unew(:, ti) = Unew(:, ti) + expeigval(:, ti - taui).*vtau;
            end
        end
        reldif = norm(Unew(:, Nt) - Ulast(:, Nt))/norm(Ulast(:, Nt));
        Ulast = Unew;
        niter = niter + 1;
    end
    if niter == Niter
        fprintf('The program has failed to reach the wanted tolerance.\n')
    end
    niter
    U(:, 2:Nt+1) = P*Unew;
end    