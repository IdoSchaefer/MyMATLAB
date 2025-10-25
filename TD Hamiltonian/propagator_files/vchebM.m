function M = vchebM(K, V, u0, dim, leftb, rightb, Ncheb)
% The function computes the vectors, which can be used for a Chebyshev
% expansion of any function of the operator, that operates on the vector u0.
% The expansion will be: f(operator)u0 = sum(a_i*M(:, i+1)), i = 0, 1, ..., Ncheb - 1 
% The coefficients a_i are f dependent, and will be computed by the
% program chebc.
% Applies for the case of a Hamiltonian composed of only x diagonal and p
% diagonal terms.
% The program will be useful when we want to compute a number of functions
% of the same operator, that operates on the same vector.
% K: a vector of the p values of the p diagonal term (kinetic energy).
% V: a vector of the x values of the x diagonal term (potential energy).
% The eigenvalue domain of the Hamiltonian is: [leftb rightb].
    M = zeros(dim, Ncheb);
    M(:, 1) = u0;
    M(:, 2) = Hpsicheb(K, V, u0, leftb, rightb);
    for k = 3:Ncheb
        M(:, k) = 2*Hpsicheb(K, V, M(:, k-1), leftb, rightb) - M(:, k-2);
    end
end