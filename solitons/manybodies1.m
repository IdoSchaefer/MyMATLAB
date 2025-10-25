function DxDv = manybodies1(t, getxv, F, r, m, N)
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
    Fvalue(N+1) = Fvalue(2);
    Dv(1) = 0;
    ibody = 2:(N-1);
    Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody))/m;
    Dv(N) = 0;
    DxDv = [Dx; Dv];
    