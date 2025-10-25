t = input('Enter the time you want to compare with:')/0.1 + 1;
pulselen = 0;
pulses = zeros(2,1);
pulsesi = 1;
i = 2;
while i <= N-1
    if allrvalues(i-1, t)<=r0
        begin = i - 1;
        while i<= N-1 && allrvalues(i, t)<=allrvalues(i-1, t)
            pulselen = pulselen + 1;
            i = i+1;
        end
        while i<= N-1 && allrvalues(i, t)>=allrvalues(i-1, t) && allrvalues(i, t)<=r0
            pulselen = pulselen + 1;
            i = i+1;
        end
        if allrvalues(i, t)>r0
            i = i + 1;
        end
        if std(allrvalues(begin:(begin + pulselen), t))>r0/200
            pulses(1, pulsesi) = begin;
            pulses(2, pulsesi) = pulselen; 
            pulsesi = pulsesi + 1;
        end
        pulselen = 0;
    else
        i=i+1;
    end
end
dim = size(pulses);
Npulse = dim(2);
npulse = 1:Npulse;
maxpulselen = max(pulses(2,:));
pulsesshp = zeros(Npulse, maxpulselen);
for ipsh = 1:Npulse
    pulsep = pulses(1, ipsh);
    pulselen = pulses(2, ipsh);
    pulsesshp(ipsh, 1:pulselen) = (allrvalues(pulsep:(pulsep + pulselen - 1), t))';
end
pulses
pulsesshp
maxpulsei = pulses(2, :);
pulsesnorm = zeros(1, Npulse);
for i = 1:Npulse
    pulsesnorm(i) = norm(pulsesshp(i,:));
end
for i = 1:(Npulse - 1)
    for j=1:(Npulse - i)
        if pulsesnorm(j)>pulsesnorm(j+1)
            ezer = pulsesnorm(j);
            pulsesnorm(j) = pulsesnorm(j+1);
            pulsesnorm(j+1) = ezer;
            ezer = npulse(j);
            npulse(j) = npulse(j+1);
            npulse(j+1) = ezer;
        end
    end
end
t2 = input('Enter the time you want to compare with the first time:')/0.1 + 1;
maxr= max(max(allrvalues));
Rezer = (allrvalues(1:(N-1), t2))' - maxr;
pulsesshp1 = pulsesshp - maxr;
maxRi = zeros(Npulse, 1);
for i = 1:Npulse
    n = npulse(i);
    maxscalmul = 0;
    for Ri=1:(N - maxpulsei(n))
        scalmulRi = sum(pulsesshp1(n, 1:maxpulsei(n)).*Rezer((1:maxpulsei(n)) + Ri -1));
        if scalmulRi>maxscalmul
            maxscalmul = scalmulRi;
            maxRi(n) = Ri;
        end
    end
    Rezer(maxRi(n):(maxRi(n)+maxpulselen-1)) = Rezer(maxRi(n):(maxRi(n)+maxpulselen-1)) - pulsesshp1(n,:);
end
maxRi