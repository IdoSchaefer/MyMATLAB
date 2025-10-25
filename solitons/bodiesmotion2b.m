figure
axtop = max(x(2*N,:));
axbot = min(x(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2]);
if N<=20
    ballsizeb = 24;
else
    ballsizeb = (20/N)*24;
end
ballsizea = sqrt(ma/mb)*ballsizeb;
for j = 1:N
    balla(j) = line(0, 0, 'marker', '.', 'markersize', ballsizea, 'color', 'r');
    ballb(j) = line(0, 0, 'marker', '.', 'markersize', ballsizeb);
end
for timei = 1:tistep:toptimei
    for j=1:N
        set(balla(j), 'xdata', x(j, timei))
        set(ballb(j), 'xdata', x(j + N, timei))
    end
    pause(tstep)
end