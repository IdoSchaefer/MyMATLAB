timei = 1:toptimei;
clf
axis([1, N-1, minra, maxra]);
xlabel('Bond index')
ylabel('r (arbitrary units)')
nrplot = line(1:(N-1), allravalues(1:(N-1), 1));
for timei = 1:tistep:toptimei
    set(nrplot, 'ydata', allravalues(1:(N-1), timei), 'color', 'r')
    pause(tstep)
end