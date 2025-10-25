function [result, all_polydeg] = bFtaylor(z, M)
    Nz = length(z);
    term = ones(Nz, 1);
    result = ones(Nz, 1);
    is_not_converged = true(Nz, 1);
    polydeg = 1;
    all_polydeg = ones(Nz, 1);
    while max(is_not_converged)
        term(is_not_converged) = z(is_not_converged).*term(is_not_converged)/(polydeg + M);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
        all_polydeg(is_not_converged) = polydeg;
    end
end