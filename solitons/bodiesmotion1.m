figure
axtop = max(x(N,:));
axbot = min(x(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2]);
xlabel('x (arbitrary units)')
if N<=32
    ballsize = 24;
else
    ballsize = (32/N)*24;
end
if size(m, 2) == 1
    for j = 1:N
        ball(j) = line(0, 0, 'marker', '.', 'markersize', ballsize);
    end
else
    maxm = max(m);
    minm = min(m);
    for j = 1:N
        ballnsz = m(j)*ballsize;
        mred = (m(j)-minm)/(maxm - minm);
        mblue = 1 - mred;
        ball(j) = line(0, 0, 'marker', '.', 'markersize', ballnsz, 'color', [mred, 0, mblue]);
    end    
end
for timei = 1:tistep:toptimei
    for j=1:N
        set(ball(j), 'xdata', x(j, timei))
    end
    pause(tstep)
end