function [result, polydeg] = Ftest4(x, t, Nt_ts, tol)
% The function of matrix F(x, t):
    Nt = length(t);
    Nx = length(x);
    minus_ixt = -1i*x*t;
%    is_big = Nt_ts*factorial(Nt_ts)*eps./abs(minus_ixt.^(Nt_ts + 1)) < tol;
    is_big = factorial(Nt_ts)*eps./abs(minus_ixt.^(Nt_ts)) < tol;
%    is_big = factorial(Nt_ts)*eps./(abs(minus_ixt.^(Nt_ts)).*abs(exp(minus_ixt))) < tol;
    result = ones(Nx, Nt);
    result(is_big) = exp(minus_ixt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_ixt(is_big);
    end
    is_not_converged = ~is_big;
    term = double(is_not_converged);
%    result(is_not_converged) = 1;
    polydeg = 1;
    while max(max(is_not_converged))
        term(is_not_converged) = minus_ixt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
    end
    result = result.*((ones(Nx, 1)*t).^Nt_ts);
end
