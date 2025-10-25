function Vnew = newguess0lb5(Vold, Nold, Nnew)
    dim = size(Vold);
    if dim(1) == 1
        Nk_free = dim(2)/2;
        Vold = Vold.';
    else
        Nk_free = dim(1)/2;
    end
    Vnew = [Vold(1:Nk_free)*sqrt(Nnew/Nold); Vold((Nk_free + 1):2*Nk_free)*(Nnew/Nold)^(1/4)];
    if dim(1) == 1
        Vnew = Vnew.';
    end    
end