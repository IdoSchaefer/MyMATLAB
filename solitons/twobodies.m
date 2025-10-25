function dxdv = twobodies(t, getxv, F, r, m)
    x = zeros(2, 1);
    v = zeros(2, 1);
    dx = zeros(2,1);
    dv = zeros(2,1);
    x(1) = getxv(1);
    v(1) = getxv(3);
    x(2) = getxv(2);
    v(2) = getxv(4);
    dx = v;
    rvalue = x(2) - x(1);
    Fvalue = subs(F, r, rvalue);
    dv(1) = -Fvalue/m;
    dv(2) = Fvalue/m;
    dxdv = [dx; dv];
    