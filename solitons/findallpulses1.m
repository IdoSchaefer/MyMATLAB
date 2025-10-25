maxr = max(max(allrvalues));
minr = min(min(allrvalues));
allpeaks = zeros(5,1, toptimei);
maxamp = max(r0 - minr, maxr - r0);
for ti = 1:toptimei
    i = 2;
    pulsesi = 1;
    while i <= N-1
        if allrvalues(i-1, ti)<=r0           
            begin = i-1;
            while i<= N-1 && allrvalues(i, ti)<=allrvalues(i-1, ti)
                i = i+1;
            end
            peak = allrvalues(i-1, ti);
            peaki = i-1;
            while i<= N-1 && allrvalues(i, ti)>=allrvalues(i-1, ti) && allrvalues(i, ti) <= r0
                i = i+1;
            end
            if (r0 - peak)/maxamp > 0.05*maxamp
                allpeaks(1, pulsesi, ti) = peak;
                allpeaks(2, pulsesi, ti) = peaki;
                allpeaks(3, pulsesi, ti) = (x(peaki, ti) + x(peaki+1, ti))/2;
                allpeaks(4, pulsesi, ti) = begin;
                allpeaks(5, pulsesi, ti) = i-1-begin;
                pulsesi = pulsesi + 1;
            end
            if i<=N-1 && allrvalues(i, ti)>r0
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end
end