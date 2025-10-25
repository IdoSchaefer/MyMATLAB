function DxDv = opencexp(t, getxv, F, Fm, r, xi, r0, m, N)
    Dx = zeros(N, 1);
    Dv = zeros(N, 1);
    ibody = 1:N;
    iforce = 2:N;
    x = getxv(ibody);
    v = getxv(N + ibody);
    xivalue = x' - (0:r0:(N-1)*r0);
    Dx(ibody) = v(ibody);
    rvalue(iforce) = x(iforce) - x(iforce-1);
    Fvalue = zeros(1, N+1);
    Fvalue(iforce) = subs(F, r, rvalue(iforce));
    Fmvalue(iforce) = subs(Fm, xi, xivalue(iforce));
    szm = size(m, 2);
    if szm == 1
        Dv(1) = (-Fvalue(2)+Fmvalue(1))/m;
        Dv(iforce) = (-Fvalue(iforce+1)+Fvalue(iforce)+Fmvalue(iforce))/m;
    else
        Dv(1) = (-Fvalue(2)+Fmvalue(1))/m(1);
        Dv(iforce) = (-Fvalue(iforce+1)+Fvalue(iforce)+Fmvalue(iforce))./m(iforce);
    end
    DxDv = [Dx; Dv];