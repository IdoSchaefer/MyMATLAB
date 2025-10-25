function [Fi, Hts] = ihchifMLSrlx(chits, miu, tsi, t_ts, Tts, Nt, H0, chimiupsi, vfilterE, ihterm, fieldw0, weight)
    firstt_ts = Nt*(-Tts) + Tts*(tsi - 1);
    dw = pi/(Nt*(-Tts));
%    wval = dw:dw:(pi/Tts - dw);
    wval = dw:dw:pi/(-Tts);
    [dim, Nt_ts] = size(chits);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
    fieldw = (1 - weight)*fieldw0 + weight*dctI(chimiupsi).*vfilterE;
    for ti = 1:Nt_ts
        tpoint = firstt_ts + t_ts(ti);
        % Performing a cosine transform of the first kind at the time point within the time step:
        fieldtpoint = sqrt(2/Nt)*(0.5*fieldw(1) + sum(fieldw(2:Nt).*cos(wval(1:(Nt-1))*tpoint)) + 0.5*fieldw(Nt+1)*cos(wval(Nt)*tpoint));
        Vtts(:, :, ti) = -miu*fieldtpoint;
    end
    Vthalf = Vtts(:, :, tmidi);
    Hts = H0 + Vthalf;
    % Calculation of the inhomogeneous fi vectors:
    for ti = 1:Nt_ts
        Fi(:, ti) = -1i*(Vtts(:, :, ti) - Vthalf)*chits(:, ti) + ihterm(:, ti, tsi);
    end    
end