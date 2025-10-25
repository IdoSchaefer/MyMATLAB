function [Fi, Hts, fieldts] = ihchifttpMLS1(chits, miu, tsi, H0, allfield, allpsi, penalM, ihterm)
    [dim, Nt_ts, Nts] = size(allpsi);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    % The indices inside the time step for allfield:
    allti_ts = (tsi - 1)*(Nt_ts - 1) + (1:Nt_ts);
    % The last index for allfield:
    allt_lasti = Nts*(Nt_ts - 1) + 1;
    for ti = 2:Nt_ts
        allfield(allti_ts(ti)) = (-imag(chits(:, ti)'*miu*allpsi(:, ti, tsi))...
            - penalM(allti_ts(ti), [1:(allti_ts(ti) - 1), (allti_ts(ti) + 1):allt_lasti])*...
            allfield([1:(allti_ts(ti) - 1), (allti_ts(ti) + 1):allt_lasti]))/penalM(allti_ts(ti), allti_ts(ti));
    end
    fieldts = allfield(allti_ts);
    for ti = 1:Nt_ts
        Vtts(:, :, ti) = -miu*fieldts(ti);
    end
    Vthalf = Vtts(:, :, tmidi);
    Hts = H0 + Vthalf;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vthalf)*chits(:, ti) + ihterm(:, ti, tsi);
    end
end