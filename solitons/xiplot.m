allxivalues = zeros(N, toptimei);
for timei = 1:toptimei
    allxivalues(:, timei) = x(1:N, timei) - (0:r0:(N-1)*r0)';
end
maxxi = max(max(allxivalues));
minxi = min(min(allxivalues));
clf
axis([1, N-1, minxi, maxxi]);
xlabel('n')
ylabel('\xi (arbitrary units)')
nxiplot = line(1:N, allxivalues(1:N, 1));
for timei = 1:tistep:toptimei
    set(nxiplot, 'ydata', allxivalues(1:N, timei))
    pause(tstep)
end