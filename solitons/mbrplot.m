timei = 1:toptimei;
allrvalues = x(2:N, timei) - x(1:(N-1), timei);
maxr = max(max(allrvalues));
minr = min(min(allrvalues));
clf
axis([1, N-1, minr, maxr]);
xlabel('Bond index')
ylabel('r (arbitrary units)')
nrplot = line(1:(N-1), allrvalues(1:(N-1), 1));
for timei = 1:tistep:toptimei
    set(nrplot, 'ydata', allrvalues(1:(N-1), timei))
    pause(tstep)
end