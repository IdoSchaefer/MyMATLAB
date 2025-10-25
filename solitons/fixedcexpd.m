function DxDv = fixedcexpd(t, getxv, F, Fm, r, xi, r0, m, N, d)
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
    Fvalue(N+1) = Fvalue(2);
    Dv(1) = 0;
    ibody = 2:(N-1);
    szm = size(m, 2);
        if szm == 1
            Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody)+Fmvalue(ibody)-d*(v(ibody))')/m;
        else
            Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody)+Fmvalue(ibody)-d*(v(ibody))')./m(ibody);
        end
    Dv(N) = 0;
    DxDv = [Dx; Dv];
    