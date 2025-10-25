function [Fi, Hts, fieldts] = ihpsittpMLS(psits, miu, tsi, H0, allfield, allchi, penalM)
    [dim, Nt_ts, Nts] = size(allchi);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    b = zeros(Nt_ts, 1);
    % The weighted penalty matrix for the new field in this propagation:
    pMnew = penalM((tsi - 1)*(Nt_ts - 1) + (1:Nt_ts), 1:((tsi - 1)*(Nt_ts - 1) + 1));
    pMnew(:, (tsi - 1)*(Nt_ts - 1) + 1) = pMnew(:, (tsi - 1)*(Nt_ts - 1) + 1)/2;
    % The weighted penalty matrix for the old field in the last propagation:
    pMold = penalM((tsi - 1)*(Nt_ts - 1) + (1:Nt_ts), (tsi*(Nt_ts - 1) + 1):(Nts*(Nt_ts - 1) + 1));
    pMold(:, 1) = pMold(:, 1)/2;
    % The weighted penalty matrix for the computation of the field in the current time step:
    pMts = penalM((tsi - 1)*(Nt_ts - 1) + (1:Nt_ts), (tsi - 1)*(Nt_ts - 1) + (1:Nt_ts));
    if tsi > 1 && tsi < Nts
        pMts(:, [1, Nt_ts]) = pMts(:, [1, Nt_ts])/2;
    elseif tsi == 1
        pMts(:, Nt_ts) = pMts(:, Nt_ts)/2;
    else
        pMts(:, 1) = pMts(:, 1)/2;
    end
    % computation of the free element b:
    for ti = 1:Nt_ts
        b(ti) = -imag(allchi(:, ti, tsi).'*miu*psits(:, ti));
    end
    b = b - pMnew*allfield(1:((tsi - 1)*(Nt_ts - 1) + 1)) - pMold*allfield((tsi*(Nt_ts - 1) + 1):(Nts*(Nt_ts - 1) + 1));
%         for tsj = 1:(tsi - 1)
%             integv = penalM((tsi - 1)*(Nt_ts - 1) + ti, (tsj - 1)*(Nt_ts - 1) + (1:Nt_ts)).*allfield((tsj - 1)*(Nt_ts - 1) + (1:Nt_ts));
%             b(ti) = b(ti) - intchebpb(integv, dt);
%         end
%         for tsj = (tsi + 1):Nts
%             integv = penalM((tsi - 1)*(Nt_ts - 1) + ti, (tsj - 1)*(Nt_ts - 1) + (1:Nt_ts)).*allfield((tsj - 1)*(Nt_ts - 1) + (1:Nt_ts));
%             b(ti) = b(ti) - intchebpb(integv, dt);
%         end
%     end
    % computation of the matrix of the coefficients:
%     Mts = penalM((tsi - 1)*(Nt_ts - 1) + (1:Nt_ts), (tsi - 1)*(Nt_ts - 1) + (1:Nt_ts));
%     for tj = 1:Nt_ts
%         Mts(:, tj) = Mts(:, tj)*vts(tj);
%     end
    fieldts = pMts\b;
%    field = -imag(allchi(:, 1, tsi)'*miu*psits(:, 1))/Epenal;
    for ti = 1:Nt_ts
        Vtts(:, :, ti) = -miu*fieldts(ti);
    end
    Vthalf = Vtts(:, :, tmidi);
    Hts = H0 + Vthalf;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vthalf)*psits(:, ti);
    end    
end