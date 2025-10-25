t = input('Enter the time:')/0.1 + 1;
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
        if i<=N-1 && allrvalues(i, t)>r0
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
pulses