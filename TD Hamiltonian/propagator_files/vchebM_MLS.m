function M = vchebM_MLS(Hcheb, u0, dim, Ncheb)
% The function computes the vectors, which can be used for a Chebyshev
% expansion of any function of the operator, that operates on the vector u0.
% The expansion will be: f(operator)u0 = sum(a_i*M(:, i+1)), i = 0, 1, ..., Ncheb - 1 
% The coefficients a_i are f dependent, and will be computed by the
% program chebc.
% The program will be useful when we want to compute a number of functions
% of the same matrix, that multiplies on the same vector.
% Hcheb: the matrix, transformed to the domain of the Chebyshev expansion:
% [-1 1]
    M = zeros(dim, Ncheb);
    M(:, 1) = u0;
    M(:, 2) = Hcheb*u0;
    for k = 3:Ncheb
        M(:, k) = 2*Hcheb*M(:, k-1) - M(:, k-2);
    end
end