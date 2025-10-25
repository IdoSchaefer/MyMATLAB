function Dpsi = TLSRKpsi(t, psit, E, chisol, Epenal)
%   chit = deval(chisol, t);
    chit = RK4interp(chisol.y(:, end:-1:1), chisol.x(end:-1:1), t);
%     Nt = length(chisol.x) - 1;
% %     if t == 0
% %         chit = chisol.y(:, Nt + 1);
% %     elseif t == chisol.x(1)
% %         chit = chisol.y(:, 1);
% %     else
%     tgi = nsmaller(chisol.x, t);
%     if chisol.x(tgi) == t
%         chit = chisol.y(:, tgi);
%     else
%         if tgi == 2
%             spi = 1:4;
%         elseif tgi == Nt+1
%             spi = (Nt-2):(Nt+1);
%         else
%             spi = (tgi-2):(tgi+1);
%         end
%         chit = NewtonIpln4(chisol.x(spi), chisol.y(:, spi), t);
%     end
    field = 1i/Epenal*(conj(chit(1))*psit(2) - conj(psit(1))*chit(2));
    Ht = [E(1)     -conj(field);
          -field   E(2)        ];
    Dpsi = -1i*Ht*psit;
end
