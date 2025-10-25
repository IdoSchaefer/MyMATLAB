function [relE, result_direct, result_taylor, allpolydeg] = bFcomparison(z, M)
% The function compares between the two ways of computation of
% \bar{f}_M(z)---the direct iterative computation, and the truncated Taylor
% computation.
% Input:
% z: A vector of the input values
% Output:
% relE: The relative error of Taylor from the direct computation.
% allpolydeg: The degree of the Taylor expansion required for convergence
% for all z values
    result_direct = bFdirect(z, M);
    [result_taylor, allpolydeg] = bFtaylor(z, M);
    relE = abs(result_taylor - result_direct)./abs(result_direct);
end