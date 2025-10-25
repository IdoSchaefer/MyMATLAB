function result = Ftest2(x, t, Nt_ts, Ncheb)
% The function of matrix F(x, t):
    minus_ixt = -1i*x*t;
    % If the argument of the exponent, -ixt, is small, the function
    % shouldn't be computed directly, because of numerical instability. We
    % use the "tail" of a Taylor expansion instead:
    if min(abs(minus_ixt)) < (Nt_ts + 1)
        term = ones(1, Ncheb);
        result = ones(1, Ncheb);
        polydeg = 1;
        minus_ixt = -1i*x*t;
        while abs(term/result) > eps
            term = minus_ixt.*term/(polydeg + Nt_ts);
            result = result + term;
            polydeg = polydeg + 1;
        end
    % If -ixt is large enough, the function is computed directly:
    else
        minus_ixt = -1i*x*t;
        result = exp(minus_ixt);
        for polyi = 1:Nt_ts
            result = polyi*(result - 1)./minus_ixt;
        end
    end
    result = result*t^Nt_ts;
end