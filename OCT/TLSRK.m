function Dpsi = TLSRK(t, psi, E, field)
%     H = diag(E);
%     H(2, 1) = -coupling(t);
%     H(1, 2) = conj(H(2, 1));
    Ht = [E(1)        -conj(field(t));
          -field(t)   E(2)            ];
    Dpsi = -1i*Ht*psi;
end
