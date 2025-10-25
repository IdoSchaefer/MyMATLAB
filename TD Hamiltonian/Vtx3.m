function Vt = Vtx3(f, x, T, Nts, Nt_ts)
% The function creates Vt for the TDH___cx solvers, with a known
% functional form of the disturbance.
% f is a function handle of the form: @(x, t). x is the grid of the x
% values.
    Vt = zeros(length(x), Nt_ts, Nts);
    szx = size(x);
    if szx(1) == 1
        x = x.';
    end
    Tts = T/Nts;
    tcheb = cos(((1:Nt_ts)*2 - 1)*pi/(2*Nt_ts));
    t_ts = 0.5*(tcheb+1)*Tts;
    for tsi = 1:Nts
        for ti = 1:Nt_ts
            Vt(:, ti, tsi) = f(x, (tsi - 1)*Tts + t_ts(ti));
        end
    end
end