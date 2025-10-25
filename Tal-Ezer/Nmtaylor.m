function Nmul = Nmtaylor(N)
    A = createA(N);
    u0 = (1:N)';
    Uexact = exact_expAv(A*5, u0);
    Nmul = 10000;
    Unum = taylor2(A, u0, 5, Nmul);
    relE = norm(Uexact - Unum)/norm(Uexact);
    while relE>1e-7
        Nmul = Nmul + 10;
        Unum = taylor2(A, u0, 5, Nmul);
        relE = norm(Uexact - Unum)/norm(Uexact);
    end
    Nmul = Nmul - 10;
    relE = 10;
    while relE>1e-7
        Nmul = Nmul + 1;
        Unum = taylor2(A, u0, 5, Nmul);
        relE = norm(Uexact - Unum)/norm(Uexact);
    end
end
    