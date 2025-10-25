function C = cheb2taylor(N)
% The function returns a matrix, contains the coefficients of the terms of the Chebychev
% polynomials, in a taylor expansion. N is the maximal degree of the
% Chebychev polynomails. The Chebychev polynomials are represented by the
% rows, and the degree of the term is represented by the columns.
    C = zeros(N+1);
    C(1, 1) = 1;
 %   C(2, 1) = 0;
    C(2, 2) = 1;
    for chebi = 3:(N+1)
        C(chebi, 1) = -C(chebi - 2, 1);
        for polyi = 2:(chebi-2)
            C(chebi, polyi) = 2*(polyi-1)*C(chebi-1, polyi-1) - C(chebi-2, polyi);
        end
        C(chebi, chebi-1) = 2*(chebi-2)*C(chebi-1, chebi-2);
        C(chebi, chebi) = 2*(chebi-1)*C(chebi-1, chebi-1);
    end
end    