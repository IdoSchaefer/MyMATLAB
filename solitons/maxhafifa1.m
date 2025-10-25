pulsesshp = [0 2:4 0;
             0 2 2 0 0];
maxscalmul = 0;
scalmulBi = zeros(2, 1);
maxBi = zeros(2, 1);
npulses = 1:2;
for Bi=1:16
    scalmulBi(1) = sum(pulsesshp(1, 1:5).*B((1:5) + Bi -1));
    Bezer = B;
    Bezer(Bi:(Bi+4)) = Bezer(Bi:(Bi+4)) - pulsesshp(1,:);
    for Bezeri=1:17
        scalmulBi(2) = sum(pulsesshp(2, 1:4).*Bezer((1:4) + Bezeri -1));
        if sum(scalmulBi)>maxscalmul
            maxscalmul = sum(scalmulBi);
            maxBi(1) = Bi;
            maxBi(2) = Bezeri;
        end
    end
end
maxBi