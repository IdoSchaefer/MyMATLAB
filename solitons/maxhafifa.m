pulsesshp = [0 2:4 0];
pulsei = 1:5;
maxscalmul = 0;
for Bi=1:16
    scalmulBi = sum(pulsesshp(pulsei).*B(pulsei + Bi -1));
    if scalmulBi>maxscalmul
        maxscalmul = scalmulBi;
        maxBi = Bi;
    end
end