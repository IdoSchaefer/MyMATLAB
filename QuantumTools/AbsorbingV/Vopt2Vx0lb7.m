function Vx = Vopt2Vx0lb7(opt_params, Nnew, Nold)
    if nargin < 3
        Nold = Nnew;
    end
    sz = size(opt_params);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        opt_params = opt_params.';
    else
        Nk_free = sz(1)/2;
    end
%     Vk = [Vopt(1:Nk_free); Vopt((Nk_free + 1):2*Nk_free).^2];
%     Vxreal = dctI0lb_2([Vk(1:Nk_free); zeros(Nnew - Nk_free, 1)])*sqrt(Nnew/Nold);
%     Vximag = [0; dstI([Vk((Nk_free + 1):2*Nk_free); zeros(Nnew - Nk_free - 1, 1)]); 0]*sqrt(Nnew/Nold) -...
%         sqrt(2/Nold)*pi/Nnew*sum(((1:Nk_free).').*Vk((Nk_free + 1):2*Nk_free))*((0:Nnew).'); 
    Vk = [opt_params(1:Nk_free); zeros(2*Nk_free, 1)];
    Vk(Nk_free + (2:2:2*Nk_free)) = opt_params((Nk_free + 1):2*Nk_free).^2;
    Vxreal = dctI0lb_2([Vk(1:Nk_free); zeros(Nnew - Nk_free, 1)]);
    Vximag = [0; dstI([Vk((Nk_free + 1):3*Nk_free); zeros(Nnew - 2*Nk_free - 1, 1)]); 0]*sqrt(Nnew/Nold) -...
        sqrt(2/Nold)*pi/Nnew*sum(((2:2:2*Nk_free).').*Vk((Nk_free + 2):2:3*Nk_free))*((0:Nnew).'); 
    Vx = Vxreal + 1i*Vximag;
    if sz(1) == 1
        Vx = Vx.';
    end
end