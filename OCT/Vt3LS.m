function Vt = Vt3LS(T, Nt)
    dt = T/Nt;
    Vt = zeros(3, 3, Nt);
    for ti = 1:Nt
        t = (ti-1)*dt;
%         Vt(:, :, ti) = [0             -exp(1i*t)      -0.1;
%                         -exp(-1i*t)   0               -exp(1i*t);
%                         -0.1          -exp(-1i*t)     0         ]*0.1;
        Vt(:, :, ti) = [0     1     0.1;
                        1     0     1;
                        0.1     1     0]*(-0.1)*cos(t);
    end
end