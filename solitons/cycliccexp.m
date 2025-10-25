function DxDv = cycliccexp(t, getxv, F, Fm, r, xi, r0, m, N)
    Dx = zeros(N, 1);
    Dv = zeros(N, 1);
    ibody = 1:N;
    ispring = 2:N;
    x = getxv(ibody);
    v = getxv(N + ibody);
    xivalue = x' - (0:r0:(N-1)*r0);
    Dx(ibody) = v(ibody);
    rvalue(1) = x(1) - (x(N) - N*r0);
    rvalue(ispring) = x(ispring) - x(ispring-1);
    Fvalue = zeros(1, N+1);
    Fvalue(ibody) = subs(F, r, rvalue(ibody));
    Fmvalue(ispring) = subs(Fm, xi, xivalue(ispring));
    ibody = 1:(N-1);
    szm = size(m, 2);
    if szm == 1    
        Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody) + Fmvalue(ibody))/m;
        Dv(N) = (-Fvalue(1) + Fvalue(N) + Fmvalue(N))/m;
    else
        Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody) + Fmvalue(ibody))./m(ibody);
        Dv(N) = (-Fvalue(1) + Fvalue(N) + Fmvalue(N))/m(N);
    end    
    DxDv = [Dx; Dv];