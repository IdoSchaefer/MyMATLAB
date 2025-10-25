function [Tv checkv] = getTv(DV, Da, Emin, Emax, NE)
    Tv = zeros(1, NE + 1);
    checkv = zeros(1, NE + 1);
    dE = (Emax - Emin)/NE;
    E = Emin;
    for Ei = 1:(NE + 1)
        [Tv(Ei) checkv(Ei)] = trans_coef(E, DV, Da);
        E = E + dE;
    end
end