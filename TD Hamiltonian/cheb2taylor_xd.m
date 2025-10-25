function C = cheb2taylor_xd(N, a, b)
% The function returns a matrix which contains the coefficients of the terms of the Chebychev
% polynomials transformed to an arbitrary domain x\in[a, b], in a taylor like expansion.
% N: The maximal degree of the Chebychev polynomials.
% C: The coefficient matrix; the Chebychev polynomials are represented by the
% rows, and the degree of the term is represented by the columns.
    D = b - a;
    S = a + b;
    C = zeros(N+1);
    C(1, 1) = 1;
    C(2, 1) = -S/D;
    C(2, 2) = 2/D;
    for chebi = 3:(N+1)
        C(chebi, 1) = -2*S/D*C(chebi - 1, 1) - C(chebi - 2, 1);
        for polyi = 2:(chebi-2)
            C(chebi, polyi) = 4/D*(polyi-1)*C(chebi-1, polyi-1) - 2*S/D*C(chebi - 1, polyi) - C(chebi-2, polyi);
        end
        C(chebi, chebi-1) = 4/D*(chebi-2)*C(chebi-1, chebi-2) - 2*S/D*C(chebi - 1, chebi - 1);
        C(chebi, chebi) = 4/D*(chebi-1)*C(chebi-1, chebi-1);
    end
end    