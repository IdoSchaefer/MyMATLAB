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
while ti<=toptimei && Npulse(ti) == maxNpulse
    ti = ti + 1;
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
    ti
    lastpeaks
    allpeaks(2, :, ti-1)
    directionb4 = sign(allpeaks(2, :, ti - 1) - lastpeaks);
    for emergedni = 1:Npulse(ti)
        distances = abs(allpeaks(3, emergedni, ti) - allpeaks(3, :, ti-1));
        [stam, mindisi] = min(distances);
        npulse(emergedni) = mindisi;
    end
    while ti<=toptimei && Npulse(ti) < maxNpulse
        ezertpeaks = allpeaks(:, :, ti);
        allpeaks(:, :, ti) = 0;
        allpeaks(:, npulse(1:Npulse(ti)), ti) = ezertpeaks(:, 1:Npulse(ti));
        ti = ti + 1;
    end
    if ti<=toptimei
        nextpeaks = allpeaks(2, :, ti);
    else
        break
    end
    for ipulses = 1:maxNpulse
        tempti = ti;
        while tempti<= toptimei && Npulse(tempti) == maxNpulse  && ...
                                            allpeaks(2, ipulses, tempti) == nextpeaks(ipulses)
            tempti = tempti + 1;
        end
        if tempti<=toptimei && Npulse(tempti) == maxNpulse
            nextpeaks(ipulses) = allpeaks(2, ipulses, tempti);
        end
    end
    directionaf = sign(nextpeaks - allpeaks(2, :, ti));
    ti
    allpeaks(2, :, ti)
    nextpeaks
    directionb4
    directionaf
    if nnz(directionb4) == maxNpulse && nnz(directionaf) == maxNpulse
        for oldni = 1:maxNpulse
            mindif = r0;
            for newni = 1:maxNpulse
                dif = abs(peaksb4(1, oldni) - allpeaks(1, newni, ti));
                if mindif>dif && directionb4(oldni) == directionaf(newni)
                    mindif = dif;
                    nmatch = newni;
                end
            end
            npulse(oldni) = nmatch;
        end
        while ti<=toptimei && Npulse(ti) == maxNpulse
            allpeaks(:, :, ti) = allpeaks(:, npulse, ti);
            ti = ti + 1;
        end
    else
        while ti<=toptimei && Npulse(ti) == maxNpulse
            ti = ti + 1;
        end
    end
end