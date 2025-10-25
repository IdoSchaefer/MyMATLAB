timei = 1:toptimei;
clf
axis([1, N-1, minrb, maxrb]);
xlabel('Bond index')
ylabel('r (arbitrary units)')
nrplot = line(1:(N-1), allrbvalues(1:(N-1), 1));
for timei = 1:tistep:toptimei
    set(nrplot, 'ydata', allrbvalues(1:(N-1), timei))
    pause(tstep)
end