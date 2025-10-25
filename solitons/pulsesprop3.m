allpeaks = zeros(5,1, toptimei);
for ti = 1:toptimei
    i = 2;
    pulsesi = 1;
    while i <= N-1
        begin = i-1;
        while i<= N-1 && allrvalues(i, ti)<=allrvalues(i-1, ti)
            i = i+1;
        end
        peak = allrvalues(i-1, ti);
        peaki = i-1;
        while i<= N-1 && allrvalues(i, ti)>=allrvalues(i-1, ti)
            i = i+1;
        end
        if maxr - peak > 0.99*(maxr - r0)
            allpeaks(1, pulsesi, ti) = peak;
            allpeaks(2, pulsesi, ti) = peaki;
            allpeaks(3, pulsesi, ti) = (x(peaki, ti) + x(peaki+1, ti))/2;
            allpeaks(4, pulsesi, ti) = begin;
            allpeaks(5, pulsesi, ti) = i-1-begin;
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
mindisi = zeros(1, 2);
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
        distances = abs(allpeaks(2, emergedni, ti) - allpeaks(2, :, ti-1));
        differences = abs(allpeaks(3, emergedni, ti) - allpeaks(3, :, ti - 1));
        [stam, mindisi(1)] = min(distances);
        distances(mindisi(1)) = N;
        [stam, mindisi(2)] = min(distances); 
        if min(differences) == differences(mindisi(1))
            npulse(emergedni) = mindisi(1);
        else
            npulse(emergedni) = mindisi(2);
        end
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
    for oldni = 1:maxNpulse
        mindif = r0;
%        differences = abs(allpeaks(1, oldni, ti) - peaksb4(1, :));
%        [stam, mind] = min(differences);
        for newni = 1:maxNpulse
            dif = abs(peaksb4(1, oldni) - allpeaks(1, newni, ti));
            if mindif>dif  && (directionb4(oldni) == directionaf(newni) || ...
                                directionaf(newni) == 0 || directionb4(oldni) == 0)                
                mindif = dif;
                nmatch = newni;
            end
        end
        npulse(oldni) = nmatch;
    end      
%        distances = abs(allpeaks(2, oldni, ti) - allpeaks(2, :, ti-1));
%        [stam, mindis(1)] = min(distances);
%        distances(mindisi(1)) = N;
%        [stam, mindisi(2)] = min(distances);
%            if 
%                npulse(oldni) = mindisi(1);
%            else
%                npulse(oldni) = mindisi(2);
%            end
    while ti<=toptimei && Npulse(ti) == maxNpulse
        allpeaks(:, :, ti) = allpeaks(:, npulse, ti);
        ti = ti + 1;
    end
%    while ti<=toptimei && Npulse(ti) == maxNpulse
%        ti = ti + 1;
%    end
end