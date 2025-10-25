function C = chebcbM(fM)
% The function computes the Chebychev coefficients of a
% function sampled at the Chebychev points that include the boundary of the
% domain.
% fM is a matrix that contains the sampled values of several functions in
% its columns.
% The Chebyshev coefficiets of each function are the corresponding columns
% of C.
    N = size(fM, 1) - 1;
    C = sqrt(2/N)*dctIM(fM);
    C([1, N+1], :) = 0.5*C([1, N+1], :);
end