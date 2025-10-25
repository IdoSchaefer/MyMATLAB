function C = cheb_rec_fun(prevCs, ~)
% The function computes the power coefficients of a particular Chebyshev
% polynomial from the power coefficients of the two previous orders.
% Input:
% prevCs: A matrix containing the power coefficients of the two previous
% orders; each order is represented by a row, and each power is represented
% by a column. Column j represents the (j-1)'th power coefficients.
% Output: A row vector of the power coefficients of the new Chebyshev
% polynomial.
    N = size(prevCs, 2);
    C = zeros(1, N);
    C(2:N) = 2*prevCs(2, 1:(N - 1));
    C(1:(N - 2)) = C(1:(N - 2)) - prevCs(1, 1:(N - 2));
end