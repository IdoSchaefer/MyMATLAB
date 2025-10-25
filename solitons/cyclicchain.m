function DxDv = cyclicchain(t, getxv, F, r, r0, m, N)
    Dx = zeros(N, 1);
    Dv = zeros(N, 1);
    ibody = 1:N;
    ispring = 2:N;
    x = getxv(ibody);
    v = getxv(N + ibody);
    Dx(ibody) = v(ibody);
    rvalue(1) = x(1) - (x(N) - N*r0);
    rvalue(ispring) = x(ispring) - x(ispring-1);
    Fvalue = zeros(1, N+1);
    Fvalue(ibody) = subs(F, r, rvalue(ibody));
    ibody = 1:(N-1);
    szm = size(m, 2);
    if szm == 1    
        Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody))/m;
        Dv(N) = (-Fvalue(1) + Fvalue(N))/m;
    else
        Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody))./m(ibody);
        Dv(N) = (-Fvalue(1) + Fvalue(N))/m(N);
    end    
    DxDv = [Dx; Dv];