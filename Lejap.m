function pout = Lejap(oldp, newp, Nout, capacity)
% The function returns the Nout Leja points, chosen from the set of points:
% newp. The old point set is oldp.
% Note: This code is an old one. There's a problem of memory management
% here. A newer version exsists, Lejap1.m (now in the folder
% NumericalMethods\interpolation).
    Nold = length(oldp);
    nnew = length(newp);
    pout = zeros(Nout, 1);
    outi = 1;
% The case of an empty oldp requires a special treatment, for the 1'st
% point - we have to choose it smartly:
    if Nold == 0        
        [stam, maxi] = max(abs(newp));
        pout(outi) = newp(maxi);
        newp = newp([1:(maxi-1), (maxi+1):nnew]);
        nnew = nnew - 1;
        Nold = 1;
        oldp = pout(outi);
        outi = outi + 1;
    end
    difmul = ones(nnew, 1);
    for oldi = 1:Nold
        difmul = difmul.*abs(newp - oldp(oldi))/capacity;
    end
    for outi = outi:Nout
        [stam, maxi] = max(difmul);
        pout(outi) = newp(maxi);
        newp = newp([1:(maxi-1), (maxi+1):nnew]);
        difmul = difmul([1:(maxi-1), (maxi+1):nnew]);
        nnew = nnew - 1;
        difmul = difmul.*abs(newp - pout(outi))/capacity;
    end
end    