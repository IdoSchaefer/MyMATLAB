function Dpsi = TILMpsi(t, psit, E, chi, Epenal)
%    field = 1i/Epenal*(conj(chi(1))*psit(2) - conj(psit(1))*chi(2));
    field = 1/Epenal*imag(conj(chi(2))*psit(1) + conj(chi(1))*psit(2));
    Ht = [E(1)     -conj(field);
          -field   E(2)        ];
    Dpsi = -1i*Ht*psit;
end
