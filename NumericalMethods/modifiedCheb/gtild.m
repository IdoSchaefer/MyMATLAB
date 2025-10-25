function result = gtild(y, alpha, beta)
% The \tilde g(y; alpha, beta) function for the nonsymmetric transformation. 
    result = asin((2*alpha*beta*y + alpha - beta)/(alpha + beta));
end