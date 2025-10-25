function [result, all_polydeg] = Ftaylor(z, t, Nt_ts)
    minus_izt = -1i*z*t;
    Nx = length(z);
    Nt = length(t);
    term = ones(Nx, Nt);
    result = ones(Nx, Nt);
    is_not_converged = true(Nx, Nt);
    polydeg = 1;
    all_polydeg = ones(Nx, Nt);
    while max(max(is_not_converged))
        term(is_not_converged) = minus_izt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
%         if polydeg>200
%             keyboard;
%         end
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
        all_polydeg(is_not_converged) = polydeg;
    end
    result = result.*((ones(Nx, 1)*t).^Nt_ts);
end