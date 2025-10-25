allvvalues = x(N+1:2*N, :);
allKvalues = 0.5*m*allvvalues.^2;
Ur0 = subs(V, r, r0);
allUvalues = subs(V, r, allrvalues) - Ur0;
pulsesK = zeros(maxNpulse, toptimei);
pulsesU = zeros(maxNpulse, toptimei);
for ti = 1:toptimei
    for pulsei = 1:maxNpulse
        pulsebegin = allpeaks(4, pulsei, ti);
        if pulsebegin ~= 0
            pulselen = allpeaks(5, pulsei, ti);
            for bodyi = pulsebegin:pulsebegin+pulselen+1
                pulsesK(pulsei, ti) = pulsesK(pulsei, ti) + allKvalues(bodyi, ti);
            end
            for ri = pulsebegin:pulsebegin+pulselen
                pulsesU(pulsei, ti) = pulsesU(pulsei, ti) + allUvalues(ri, ti);
            end    
        end
    end
end
pulsesE = pulsesK + pulsesU;