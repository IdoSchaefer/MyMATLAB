function P = polapp(x, a, y)
    P = 0;
    N = length(x);
    for i = 1:N
        term = a(i);
        for j = 1:i-1
            term = term*(y - x(j));
        end
        P = P + term;
    end
end