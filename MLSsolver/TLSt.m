function Dc = TLSt(t, c, E, coupling)
    H = diag(E);
    H(1, 2) = coupling(t);
    H(2, 1) = conj(H(1, 2));
    Dc = -1i*H*c;
end
