function [Fi, Hts] = ihchifMLS(chits, miu, tsi, t_ts, Tts, Nt, H0, chimiupsi, vfilterE, ihterm)
    firstt_ts = Nt*(-Tts) + Tts*(tsi - 1);
    dw = pi/(Nt*(-Tts));
%    wval = dw:dw:(pi/Tts - dw);
    wval = dw:dw:pi/(-Tts);
    [dim, Nt_ts] = size(chits);
    Vtts = zeros(dim, dim, Nt_ts);
    Fi = zeros(dim, Nt_ts);
    tmidi = round(Nt_ts/2);
%    fieldw = dct(chimiupsi(1, :)).*vfilterE;
    fieldw = dctI(chimiupsi).*vfilterE;
    for ti = 1:Nt_ts
%         fieldtpoint = sqrt(1/Nt)*fieldw(1) + sqrt(2/Nt)*sum(fieldw(2:Nt).*cos(wval*(firstt_ts + t_ts(ti) + Tts/2)));
% % The argument of cos is simply wt. t is chosen in accordance to the time
% % grid of the dct, which starts in dt/2, and ends in T-dt/2.
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