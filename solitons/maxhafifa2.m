pulsesshp = [2:4;
             2 2 0;
             5 5 7;
             4:-1:2];
dim = size(pulsesshp);
maxpulselen = dim(2);
Npulse = dim(1);
npulse = 1:Npulse;
maxpulsei = zeros(1, Npulse);
for i = 1:Npulse
    maxpulsei(i) = nnz(pulsesshp(i, :));
end
pulsesnorm = zeros(1, Npulse);
for i = 1:Npulse
    pulsesnorm(i) = norm(pulsesshp(i,:));
end
for i = 1:(Npulse - 1)
    for j=1:(Npulse - i)
        if pulsesnorm(j)<pulsesnorm(j+1)
            ezer = pulsesnorm(j);
            pulsesnorm(j) = pulsesnorm(j+1);
            pulsesnorm(j+1) = ezer;
            ezer = npulse(j);
            npulse(j) = npulse(j+1);
            npulse(j+1) = ezer;
        end
    end
end
Bezer = B;
maxBi = zeros(Npulse, 1);
for i = 1:Npulse
    n = npulse(i);
    maxscalmul = 0;
    for Bi=1:16
        scalmulBi = sum(pulsesshp(n, 1:maxpulsei(n)).*Bezer((1:maxpulsei(n)) + Bi -1));
        if scalmulBi>maxscalmul
            maxscalmul = scalmulBi;
            maxBi(n) = Bi;
        end
    end
    Bezer(maxBi(n):(maxBi(n)+maxpulselen-1)) = Bezer(maxBi(n):(maxBi(n)+maxpulselen-1)) - pulsesshp(n,:);
end
maxBi