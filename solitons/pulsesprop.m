allpeaks = zeros(3,1, toptimei);
for ti = 1:toptimei
    i = 2;
    pulsesi = 1;
    while i <= N-1
        while i<= N-1 && allrvalues(i, ti)<=allrvalues(i-1, ti)
            i = i+1;
        end
        peak = allrvalues(i-1, ti);
        peaki = i-1;
        while i<= N-1 && allrvalues(i, ti)>=allrvalues(i-1, ti)
            i = i+1;
        end
        if (r0 - peak)/A > 0.1
            allpeaks(1, pulsesi, ti) = peak;
            allpeaks(2, pulsesi, ti) = peaki;
            allpeaks(3, pulsesi, ti) = (x(peaki, ti) + x(peaki+1, ti))/2;
            pulsesi = pulsesi + 1;
        end
    end
end
Npulse = zeros(1, toptimei);
for ti=1:toptimei
    Npulse(ti) = nnz(allpeaks(1, :, ti));
end
maxNpulse = max(Npulse);
npulse = zeros(1, maxNpulse);
ti = 1;
while Npulse(ti) < maxNpulse
    ti = ti + 1;
end    
tmaxNpulse = ti;
%lastpeaks = zeros(1, 4);
while ti<=toptimei && Npulse(ti) == maxNpulse
    ti = ti + 1;
%    if allpeaks(2, :, ti - 1) ~= allpeaks(2, :, ti-2)
%        for pulsesi = 1:maxNpulse
%            if allpeaks(2, pulsesi, ti-1) ~= allpeaks(2, pulsesi, ti-2)
%                lastpeaks(pulsesi) = allpeaks(2, pulsesi, ti-2);
%            end
%        end   
%    end
end
while ti <= toptimei    
    peaksb4 = allpeaks(: ,: ,ti - 1);
    lastpeaks = allpeaks(2, :, ti - 1);
    for ipulses = 1:maxNpulse
        tempti = ti - 1;
        while tempti >= 1 && Npulse(tempti) == maxNpulse  && ...
                                    allpeaks(2, ipulses, tempti) == lastpeaks(ipulses)
            tempti = tempti - 1;
        end
        if Npulse(tempti) == maxNpulse
            lastpeaks(ipulses) = allpeaks(2, ipulses, tempti);
        end
    end
    directionb4 = sign(allpeaks(2, :, ti - 1) - lastpeaks);
%        for ipulses = 1:maxNpulse
%            if directionb4(ipulses) == 0
%                mindistance = N;
%                imindistance = 
%                for icompulse = 1:maxNpulse
%                    if ipulses ~= icompulse && 
%                                        abs(peaksb4(2, ipulses) - peaksb4(2, icompulse)) < mindistance
%                        directionb4(ipulses) == directionb4(                        
    wallcol = zeros(1, maxNpulse - 1);
%check direction change after collision with the boundary.
    ti
    allpeaks(2, :, ti - 1)
    while ti<=toptimei && Npulse(ti) < maxNpulse
        for ipulses=1:(maxNpulse-1)
            if (allpeaks(2, ipulses, ti) >= 1 && allpeaks(2, ipulses, ti) <= 5) ...
                        || (allpeaks(2, ipulses, ti) <= N-1 && allpeaks(2, ipulses, ti) >= N-6)
                wallcol(ipulses) = 1;
            end
        end
        ti = ti + 1;
    end
%assumption (wrong): the number of pulses doesn't change during the period of
%combined pulses.
    if nnz(wallcol) > 0
        lastcompeaks = allpeaks(2, :, ti - 1);
        Nlast = nnz(lastcompeaks);
        separate = zeros(2, Nlast);
        sepdis = ones(2, Nlast)*N*r0;
        for ilast = 1:Nlast
            if wallcol(ilast)
                for ipulses=1:maxNpulse
                    if sepdis(1, ilast) > abs(lastcompeaks(ilast) - allpeaks(2, ipulses, ti))
                        sepdis(2, ilast) = sepdis(1, ilast);
                        separate(2, ilast) = separate(1, ilast);
                        sepdis(1, ilast) = abs(lastcompeaks(ilast) - allpeaks(2, ipulses, ti));
                        separate(1, ilast) = ipulses;
                    elseif sepdis(2, ilast) > abs(lastcompeaks(ilast) - allpeaks(2, ipulses, ti))
                        sepdis(2, ilast) = abs(lastcompeaks(ilast) - allpeaks(2, ipulses, ti));
                        separate(2, ilast) = ipulses;                        
                    end
                end
            end
        end
    end
    nextpeaks = allpeaks(2, :, ti);
    for ipulses = 1:maxNpulse
        tempti = ti;
        while tempti<= toptimei && Npulse(tempti) == maxNpulse  && ...
                                            allpeaks(2, ipulses, tempti) == nextpeaks(ipulses)
            tempti = tempti + 1;
        end
        if Npulse(tempti) == maxNpulse
            nextpeaks(ipulses) = allpeaks(2, ipulses, tempti);
        end
    end
    directionaf = sign(nextpeaks - allpeaks(2, :, ti));
    if nnz(wallcol) > 0
        for i=1:2*Nlast
            if separate(i)
                directionaf(separate(i)) = -directionaf(separate(i));
            end
        end
    end
    ti
    allpeaks(2, :, ti)
    lastpeaks
    nextpeaks
    directionb4
    directionaf
    wallcol
    if nnz(directionb4) == maxNpulse && nnz(directionaf) == maxNpulse
        for oldni = 1:maxNpulse
            mindif = r0;
            for newni = 1:maxNpulse
                dif = abs(peaksb4(1, oldni) - allpeaks(1, newni, ti));
                if mindif>dif %&& directionb4(oldni) == directionaf(newni)
                    mindif = dif;
                    nmatch = newni;
                elseif mindif == dif && directionb4(oldni) == directionaf(newni)
                    nmatch = newni;    
                end
            end
            npulse(oldni) = nmatch;
        end
        while ti<=toptimei && Npulse(ti) == maxNpulse
            allpeaks(:, :, ti) = allpeaks(:, npulse, ti);
            ti = ti + 1;
%        if allpeaks(2, :, ti - 1) ~= allpeaks(2, :, ti-2)
%            for pulsesi = 1:maxNpulse
%                if allpeaks(2, pulsesi, ti-1) ~= allpeaks(2, pulsesi, ti-2)
%                    lastpeaks(pulsesi) = allpeaks(2, pulsesi, ti-2);
%                end
%            end   
%        end
        end
    else
        while ti<=toptimei && Npulse(ti) == maxNpulse
            ti = ti + 1;
        end
    end
end