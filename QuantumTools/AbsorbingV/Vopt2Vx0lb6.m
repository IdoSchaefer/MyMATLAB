function Vx = Vopt2Vx0lb6(Vopt, Nnew, Nold)
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
    Vximag = [0; dstI([Vk((Nk_free + 1):2*Nk_free); zeros(Nnew - Nk_free - 1, 1)]); 0]*sqrt(Nnew/Nold) -...
        sqrt(2/Nold)*pi/Nnew*sum(((1:Nk_free).').*Vk((Nk_free + 1):2*Nk_free))*((0:Nnew).'); 
    Vx = Vxreal + 1i*Vximag;
    if sz(1) == 1
        Vx = Vx.';
    end
end