t = input('Enter the time:')/0.1 + 1;
pulselen = 0;
pulses = zeros(2,1);
pulsesi = 1;
i = 2;
while i <= N-1
    begin = i-1;
    while i<= N-1 && allrvalues(i, t)<=allrvalues(i-1, t)
        pulselen = pulselen + 1;
        i = i+1;
    end
    while i<= N-1 && allrvalues(i, t)>=allrvalues(i-1, t)
        pulselen = pulselen + 1;
        i = i+1;
    end
    if std(allrvalues(begin:(begin + pulselen), t))>r0/200
        pulses(1, pulsesi) = begin;
        pulses(2, pulsesi) = pulselen; 
        pulsesi = pulsesi + 1;
    end
    pulselen = 0;
end
dim = size(pulses);
Npulse = dim(2);
maxpulseslen = max(pulses(2,:));
pulsesshp = zeros(Npulse, maxpulseslen);
for ipsh = 1:Npulse
    pulsep = pulses(1, ipsh);
    pulselen = pulses(2, ipsh);
    pulsesshp(ipsh, 1:pulselen) = (allrvalues(pulsep:(pulsep + pulselen - 1), t))';
end
pulses
pulsesshp