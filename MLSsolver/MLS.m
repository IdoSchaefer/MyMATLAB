function Dc = MLS(t, c, E, coupling)
% The fuction assumes the same coupling between any two close eigenstates, in the order 
% of the vector E
    H = diag(E);
    N = length(E);
    for leveli = 1:N-1
        H(leveli, leveli+1) = coupling(t);
        H(leveli+1, leveli) = conj(H(leveli, leveli+1));
    end
    Dc = -1i*H*c;
end
