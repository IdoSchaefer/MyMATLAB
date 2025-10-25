function [result, relE] = test_c2t(f, max_x, x, N)
    Ctaylor = chexp2tay(f, max_x, N)
    result = 0;
    for polyi = 1:N
        result = result + Ctaylor(polyi)*x^(polyi-1)/factorial(polyi-1);
    end
    relE = (result - f(x))/f(x);
end
    