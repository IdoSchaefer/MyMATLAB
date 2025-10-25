function Dchi = TLSTDchif(t, chit, E, psisol, Epenal, target, twheight)
% Intended for a target of an analytical functional form.
% target is a function handle.
%    psit = deval(psisol, t);
    psit = RK4interp(psisol.y, psisol.x, t);
    field = 1i/Epenal*(conj(chit(1))*psit(2) - conj(psit(1))*chit(2));
    Ht = [E(1)     -conj(field);
          -field   E(2)        ];
    tart = target(t);
    Dchi = -1i*Ht*chit - twheight(t)*(tart'*psit)*tart;
end
