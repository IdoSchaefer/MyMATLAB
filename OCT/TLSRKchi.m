function Dchi = TLSRKchi(t, chit, E, psisol, Epenal)
%    psit = deval(psisol, t);
    psit = RK4interp(psisol.y, psisol.x, t);
    field = 1i/Epenal*(conj(chit(1))*psit(2) - conj(psit(1))*chit(2));
    Ht = [E(1)     -conj(field);
          -field   E(2)        ];
    Dchi = -1i*Ht*chit;
end
