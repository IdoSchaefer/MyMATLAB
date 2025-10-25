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
        if (r0 - peak)/A > 0.05*A
            allpeaks(1, pulsesi, ti) = peak;
            allpeaks(2, pulsesi, ti) = peaki;
            allpeaks(3, pulsesi, ti) = (x(peaki, ti) + x(peaki+1, ti))/2;
            allpeaks(4, pulsesi, ti) = begin;
            allpeaks(5, pulsesi, ti) = i-1-begin;
            pulsesi = pulsesi + 1;
        end
    end
end