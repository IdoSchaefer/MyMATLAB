function result = sqabs(A)
% The function returns the square of the absolute value of the terms of
% a complex array A.
    result = conj(A).*A;
end