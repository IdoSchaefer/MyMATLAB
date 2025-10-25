function Vt = coupM4(coupling, T, Nts, Nt_ts)
% The function creates Vt for the TDH___c solvers, for a TLS, with a known
% functional form of the disturbance, represented by the function handle coupling.
    Vt = zeros(2, 2, Nts, Nt_ts);
    Tts = T/Nts;
    tcheb = cos(((1:Nt_ts)*2 - 1)*pi/(2*Nt_ts));
    t_ts = 0.5*(tcheb+1)*Tts;
    for tsi = 1:Nts
        for ti = 1:Nt_ts
            Vt(1, 2, tsi, ti) = coupling((tsi - 1)*Tts + t_ts(ti));
            Vt(2, 1, tsi, ti) = conj(Vt(1, 2, tsi, ti));
        end
    end
end