function result = exvalw(A, psi, w)
% The function returns the expectation value of the operator A, for an unnormalized wave
% function, psi, in a grid with a weight function vector w.
% A is a matrix, and psi is a vector.
    result = (psi'*A*(w.*psi))/(psi'*(w.*psi));
end