function result = Ftest1(x, t, Nt_ts, Nkr)
% The function of matrix F(x, t):
    minus_ixt = -1i*x*t;
    % If the argument of the exponent, -ixt, is small, the function
    % shouldn't be computed directly, because of numerical instability. We
    % use the "tail" of a Taylor expansion instead:
    is_small = abs(minus_ixt) < Nt_ts + 1;
    is_big = ~is_small;
    term = double(is_small);
    result = double(is_small);
    polydeg = 1;
    while max(max(is_small))
        term(is_small) = minus_ixt(is_small).*term(is_small)/(polydeg + Nt_ts);
        result(is_small) = result(is_small) + term(is_small);
        polydeg = polydeg + 1;
        is_small(is_small) = abs(term(is_small)./result(is_small)) > eps;
    end
    % If -ixt is large enough, the function is computed directly:
    result(is_big) = exp(minus_ixt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_ixt(is_big);
    end
    result = result.*((ones(Nkr, 1)*t).^Nt_ts);
end
