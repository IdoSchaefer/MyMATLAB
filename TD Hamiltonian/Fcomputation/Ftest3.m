function [result, polydeg] = Ftest3(x, t, Nt_ts, tol)
% The function of matrix F(x, t):
    Nt = length(t);
    Nx = length(x);
    minus_ixt = -1i*x*t;
    is_big = true(Nx, Nt);
    % First, we try computing the function directly:
    result = exp(minus_ixt);
    result_minus1 = zeros(Nx, Nt);
    %for polyi = 1:Nt_ts
    polyi = 1;
    while polyi<=Nt_ts && max(max(is_big))
        result_minus1(is_big) = result(is_big) - 1;
        % A check if -ixt is large enough to be computed directly without roundoff errors:
        is_big(is_big) = eps./abs(result_minus1(is_big))<tol;
        result(is_big) = polyi*(result_minus1(is_big))./minus_ixt(is_big);
        polyi = polyi + 1;
    end
    % If the argument of the exponent, -ixt, is small, the function
    % shouldn't be computed directly, because of roundoff errors. We
    % use the "tail" of a Taylor expansion instead:
    is_not_converged = ~is_big;
    term = double(is_not_converged);
    result(is_not_converged) = 1;
    polydeg = 1;
    while max(max(is_not_converged))
        term(is_not_converged) = minus_ixt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
    end
    result = result.*((ones(Nx, 1)*t).^Nt_ts);
end
