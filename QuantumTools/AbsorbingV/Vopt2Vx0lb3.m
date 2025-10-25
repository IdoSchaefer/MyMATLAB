function Vx = Vopt2Vx0lb3(Vopt, Nnew)
    sz = size(Vopt);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        Vopt = Vopt.';
    else
        Nk_free = sz(1)/2;
    end
    Vxreal = dctI([Vopt(1:Nk_free); -(Vopt(1)/2 + sum(Vopt(2:Nk_free))); zeros(Nnew - Nk_free, 1)]);
    dctVkimag = dctI([Vopt((Nk_free + 1):2*Nk_free); -(Vopt(Nk_free + 1)/2 + sum(Vopt((Nk_free + 2):2*Nk_free))); zeros(Nnew - Nk_free, 1)]); 
    Vximag = -dctVkimag.^2;
    Vx = Vxreal + 1i*Vximag;
    if sz(1) == 1
        Vx = Vx.';
    end
end