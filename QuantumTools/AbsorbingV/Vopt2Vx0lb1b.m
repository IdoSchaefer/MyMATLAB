function Vx = Vopt2Vx0lb1b(Vopt, Nnew)
    sz = size(Vopt);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        Vopt = Vopt.';
    else
        Nk_free = sz(1)/2;
    end
    Vxreal = dctI0lb_2([Vopt(1:Nk_free); zeros(Nnew - Nk_free, 1)]);
    Vximag = dctI0lb_2([Vopt((Nk_free + 1):2*Nk_free); zeros(Nnew - Nk_free, 1)]); 
    Vx = Vxreal + 1i*Vximag;
    if sz(1) == 1
        Vx = Vx.';
    end
end