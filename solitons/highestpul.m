highestpeak = zeros(5, 1, toptimei);
allpeaks(allpeaks==0) = NaN;
for ti = 1:toptimei
    [stam, nhighest] = min(allpeaks(1, :, ti));
    highestpeak(:, 1, ti) = allpeaks(:, nhighest, ti);
end
highestpeak(isnan(highestpeak)) = 0;
allpeaks = highestpeak;
maxNpulse = 1;
tmaxNpulse = 1;