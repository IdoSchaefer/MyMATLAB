figure
axtop = max(x(N,:));
axbot = min(x(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2]);
if N<=32
    ballsize = 24;
else
    ballsize = (32/N)*24;
end
for j = 1:N
    ball(j) = line(0, 0, 'marker', '.', 'markersize', ballsize);
end
for timei = 1:toptimei
    for j=1:N
        set(ball(j), 'xdata', x(j, timei))
    end
    pause(0.1)
end