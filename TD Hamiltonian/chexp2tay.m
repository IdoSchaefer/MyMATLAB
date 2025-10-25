function Ctaylor = chexp2tay(f, max_x, N)
    Ctaylor = zeros(1, N);
    Ccheb = chebc(f, 0, max_x, N);
    Cc2t = cheb2taylor(N-1);
    Ctaylor(N) = (2/max_x)^(N-1)*Cc2t(N, N)*Ccheb(N);
    for polyi = (N-1):-1:1
        Ctaylor(polyi) = Cc2t(polyi, polyi)*Ccheb(polyi);
        for polyj = (polyi+1):N
            Ctaylor(polyi) = Ctaylor(polyi) + Cc2t(polyj, polyi)*Ccheb(polyj) -...
            (max_x/2)^(polyj-1)/factorial(polyj - polyi)*Ctaylor(polyj);
        end
        Ctaylor(polyi) = (2/max_x)^(polyi-1)*Ctaylor(polyi);
    end