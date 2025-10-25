function U = TDHcheb_tsf(H0, minE0, maxE0, Vt, ui, T, Nt, Nt_ts, Ncheb, tol)
% Vt is a function of the form: @(x, t).
% Intended for the case when Vt is a function of the variable of
% ui, and has a known analytical form.
% x is the grid of ui.
    dt = T/Nt;
    Tts = Nt_ts*dt;
    Nts = Nt/Nt_ts;
    dim = length(ui);
    I = eye(dim);
    U = zeros(dim, Nt + 1);
    U(:, 1) = ui;
    Ccheb = zeros(Ncheb, Nt_ts);
    for ti = 1:Nt_ts
        ft = @(x) exp(-1i*x*ti*dt);
        Ccheb(:, ti) = chebc(ft, minE0, maxE0, Ncheb).';
    end
    H0cheb = (2*H0 - (minE0 + maxE0)*I)/(maxE0 - minE0);
    Niter = 10000;
    for tsi = 1:Nts        
% Computation of U0 - the solution within the time step, for H = H0:
        v1 = U(:, (tsi-1)*Nt_ts + 1);
        v2 = H0cheb*v1;
        U0 = v1*Ccheb(1,:) + v2*Ccheb(2, :);
        for k = 3:Ncheb
            vk = 2*H0cheb*v2 - v1;
            U0 = U0 + vk*Ccheb(k, :);
            v1 = v2;
            v2 = vk;
        end
        Ulast = U0; 
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            Unew = U0;
            for taui = 1:Nt_ts
                v1 = -1i*Vt(x, (tsi - 1)*Tts + taui*dt).*Ulast(:, taui)*dt;
                v2 = H0cheb*v1;
                % Adding the term of ti = taui to the sum, which is simply v1:
                Unew(:, taui) = Unew(:, taui) + v1;
%             explanation: ti = (taui+1):Nt_ts
%                                             ti - taui = 1:(Nt_ts-taui)
                Unew(:, (taui+1):Nt_ts) = Unew(:, (taui+1):Nt_ts) + v1*Ccheb(1, 1:(Nt_ts-taui)) + v2*Ccheb(2, 1:(Nt_ts-taui));
                for k = 3:Ncheb
                    vk = 2*H0cheb*v2 - v1;
                    Unew(:, (taui+1):Nt_ts) = Unew(:, (taui+1):Nt_ts) + vk*Ccheb(k, 1:(Nt_ts-taui));
                    v1 = v2;
                    v2 = vk;
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
end