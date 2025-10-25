Ei = 1;
newpE = zeros(1, 9200/200 + 1);
for ti=1:200:9201
    newpE(Ei) = mean(pulsesE(ti:ti+49));
    Ei = Ei + 1;
end
DE=(newpE(2:9200/200+1) - newpE(1:9200/200))/2;