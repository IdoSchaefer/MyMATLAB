function Vx = Vopt2Vx0lb5(Vopt, Nnew, Nold)
    if nargin < 3
        Nold = Nnew;
    end
    sz = size(Vopt);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        Vopt = Vopt.';
    else
        Nk_free = sz(1)/2;
    end
    Vk = [Vopt(1:Nk_free); Vopt((Nk_free + 1):2*Nk_free).^2];
    Vxreal = dctI0lb_2([Vk(1:Nk_free); zeros(Nnew - Nk_free, 1)])*sqrt(Nnew/Nold);
    Vximag = dctI0lb_2([Vk((Nk_free + 1):2*Nk_free); zeros(Nnew - Nk_free, 1)])*sqrt(Nnew/Nold); 
    Vx = Vxreal + 1i*Vximag;
    if sz(1) == 1
        Vx = Vx.';
    end
end