function intv = interpct(v, Nint)
% The function computes the dctI interpolation of a signal, in Nint + 1
% equidistant points that include the boundaries. It is similiar to the
% interpft function of MATLAB.
% v: The original signal
% intv: The interpolated signal
    szv = size(v);
    if szv(1) == 1
        Nv = szv(2) - 1;
        v = v.';
    else
        Nv = szv(1) - 1;
    end
    dctv = dctI(v);
    intv = dctI([dctv(1:Nv); dctv(Nv + 1)/2; zeros(Nint - Nv, 1)])*sqrt(Nint/Nv);
    if szv(1) == 1
        intv = intv.';
    end
end