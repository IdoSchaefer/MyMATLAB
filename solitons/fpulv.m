topti = toptimei-tmaxNpulse + 1;
pulsesv = zeros(maxNpulse, topti);
pulsesn(1:maxNpulse, 1:topti) = allpeaks(2, :, tmaxNpulse:toptimei);
for pulsei = 1:maxNpulse
    ti = 1;
    while ti <= topti
        %to find the next change in n:
        nextnti = ti + 1;
        while nextnti<=topti && pulsesn(pulsei, ti) == pulsesn(pulsei, nextnti)
            nextnti = nextnti + 1;
        end
        if nextnti <= topti && pulsesn(pulsei, nextnti) ~= 0 && pulsesn(pulsei, ti) ~= 0
            vpulse = (pulsesn(pulsei, nextnti) - pulsesn(pulsei, ti))/((nextnti - ti)*dt);
        else
            vpulse = 0;
        end
        while ti<=topti && ti <= nextnti
            pulsesv(pulsei, ti) = vpulse;
            ti = ti + 1;
        end
    end
end     