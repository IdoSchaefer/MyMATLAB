function Dc = MLScheck(t, c, E, coupling)
% The fuction assumes the same coupling between any two close eigenstates, in the order 
% of the vector E
    H = diag(E);
%    N = length(E);
    H(1, 2) = coupling(t);
    H(2, 1) = conj(H(1, 2));
    H(3, 2) = coupling(t);
    H(2, 3) = conj(H(1, 2));
    Dc = -1i*H*c;
end
