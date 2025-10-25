function [U, matvecs] = NLHroie(H0, Vt, ui, x, T, Nts, Nt_ts, tol)
    Tts = T/Nts;
    dim = length(ui);
    U = zeros(dim, Nts + 1);
    U(:, 1) = ui;
    szx = size(x);
    if szx(1) == 1
        x = x.';
    end
    matvecs = 0;
    % Computation of the Imn:
    S = zeros(Nt_ts);
    theta = ((1:Nt_ts)*2 - 1)*pi/(2*Nt_ts);
% The rows of C represent the degree+1 of the Chebychev polynomials. The
% columns represent the different Chebychev points.
    C = cos((0:Nt_ts-1)'*theta);
% S is arranged like C.
    S(1, :) = 1 + cos(theta);
    S(2, :) = -sin(theta).^2/2;
    S(3:Nt_ts, :) = (((-ones(Nt_ts-2, 1)).^((3:Nt_ts)'))*ones(1, Nt_ts)...
        - (ones(Nt_ts-2, 1)*cos(theta)).*C(3:Nt_ts, :) - ((2:Nt_ts-1)'*sin(theta)).*sin((2:Nt_ts-1)'*theta))./...
        ((((2:Nt_ts-1).').^2)*ones(1, Nt_ts) - 1);
%    S(1, :) = 1 - cos(theta);
%    S(2, :) = sin(theta).^2/2;
%    S(3:Nt_ts, :) = ((ones(Nt_ts-2, 1)*cos(theta)).*C(3:Nt_ts, :) - 1 + ((2:Nt_ts-1)'*sin(theta)).*sin((2:Nt_ts-1)'*theta))./...
%        ((((2:Nt_ts-1).').^2)*ones(1, Nt_ts) - 1);
    I = (C.'*S - C(1, :).'*S(1, :)/2)/Nt_ts;
    tcheb = cos(theta);
    t_ts = Tts/2*(1+tcheb);
    Ulast = zeros(dim, Nt_ts);
    Unew = zeros(dim, Nt_ts);
    Zeta = zeros(dim, Nt_ts);
    Niter = 1000;
    for tsi = 1:Nts
% First guess:
        for ti = 1:Nt_ts
            Ulast(:, ti) = exp(-1i*exval(H0 + diag(Vt(U(:, tsi), x, (tsi-1)*Tts + t_ts(ti))), ui)*t_ts(ti))*ui;
        end
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Computation of the Zeta vectors:
            for ti = 1:Nt_ts
                Zeta(:, ti) = -1i*Tts*(H0 + diag(Vt(Ulast(:, ti), x, (tsi-1)*Tts + t_ts(ti))))*Ulast(:, ti);
            end
            % Computation of the new wave functions:
            for tj = 1:Nt_ts
                Unew(:, tj) = U(:, tsi);
                for ti = 1:Nt_ts
                    Unew(:, tj) = Unew(:, tj) + Zeta(:, ti)*I(ti, tj);
                end
            end
            reldif = norm(Unew(:, 1) - Ulast(:, 1))/norm(Ulast(:, 1));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            fprintf('The program has failed to reach the desired tolerance.\n')
        end        
        matvecs = matvecs + niter*Nt_ts;
        % Computation of the wave function at the edge of the time step:
        Cchebt = dct(Unew.').'/sqrt(Nt_ts);
        Cchebt(:, 2:Nt_ts) = Cchebt(:, 2:Nt_ts)*sqrt(2);
        U(:, tsi+1) = sum(Cchebt, 2);
    end
end