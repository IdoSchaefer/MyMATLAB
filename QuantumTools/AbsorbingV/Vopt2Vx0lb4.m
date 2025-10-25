function Vx = Vopt2Vx0lb4(Vopt, Nnew)
    sz = size(Vopt);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        Vopt = Vopt.';
    else
        Nk_free = sz(1)/2;
    end
    if Nk_free<Nnew
        Vxreal = dctI([Vopt(1:Nk_free); -(Vopt(1)/2 + sum(Vopt(2:Nk_free))); zeros(Nnew - Nk_free, 1)]);
    else
        Vxreal = dctI([Vopt(1:Nk_free); -2*(Vopt(1)/2 + sum(Vopt(2:Nk_free)))]);
    end
    dctVoptimag = dctI([Vopt((Nk_free + 1):2*Nk_free); -(Vopt(Nk_free + 1)/2 + sum(Vopt((Nk_free + 2):2*Nk_free))); zeros(2*Nnew - Nk_free, 1)]); 
    Vximag = -dctVoptimag(1:(Nnew + 1)).^2;
    Vx = Vxreal + 1i*Vximag;
    if sz(1) == 1
        Vx = Vx.';
    end
end