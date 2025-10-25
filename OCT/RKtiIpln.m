function tiIpln = RKtiIpln(t, Nt, dt)
    midti = t/dt + 1;
    if midti>=2 && midti<Nt
        floorti = floor(midti);
        tiIpln = (floorti - 1):(floorti + 2);
    elseif midti<2
        tiIpln = 1:4;
    else
        tiIpln = (Nt - 2):(Nt + 1);
    end
end