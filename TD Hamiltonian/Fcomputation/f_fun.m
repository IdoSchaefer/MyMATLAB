function result = f_fun(z, t, Nt_ts, tol)
% The function f(z, t):
    Nt = length(t);
    Nz = length(z);
    zt = z*t;
    % Condition for estimating if f(z, t) should be computed directly or by
    % a "tail" of a Taylor expansion (see supplementary material):
    is_big = factorial(Nt_ts)*eps./abs(zt.^(Nt_ts)) < tol;
    result = ones(Nz, Nt);
    % First, we compute f(z, t)/(t^Nt_ts), which is a function of zt.
    % A direct computation for large arguments:
    result(is_big) = exp(zt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./zt(is_big);
    end
    % Computation by a Taylor form for small arguments:
    is_not_converged = ~is_big;
    term = double(is_not_converged);
    polydeg = 1;
    while max(max(is_not_converged))
        term(is_not_converged) = zt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
    end
    % Obtaining the required function f(z,t):
    result = result.*((ones(Nz, 1)*t).^Nt_ts);
end
