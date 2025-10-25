function DxDv = openchain(t, getxv, F, r, m, N)
    Dx = zeros(N, 1);
    Dv = zeros(N, 1);
    ibody = 1:N;
    iforce = 2:N;
    x = getxv(ibody);
    v = getxv(N + ibody);
    Dx(ibody) = v(ibody);
    rvalue(iforce) = x(iforce) - x(iforce-1);
    Fvalue = zeros(1, N+1);
    Fvalue(iforce) = subs(F, r, rvalue(iforce));
    szm = size(m, 2);
    if szm == 1
        Dv(1) = -Fvalue(2)/m;
        Dv(iforce) = (-Fvalue(iforce+1)+Fvalue(iforce))/m;
    else
        Dv(1) = -Fvalue(2)/m(1);
        Dv(iforce) = (-Fvalue(iforce+1)+Fvalue(iforce))./m(iforce);
    end
    DxDv = [Dx; Dv];