function Ctaylor = Cpoly2Ctaylor(Cpoly, Cp2t)
% The function computes the Taylor-like coefficients from the coefficients
% of a polynomial set, using the conversion coefficients Cp2t.
% Input:
% Cpoly: A matrix containing the vector coefficients of the polynomial set
% in separate columns.
% Cp2t: A matrix containing the conversion coefficients from the
% polynomial set to the Taylor polynomials (see, for example, r2Taylor4).
% Output:
% Ctaylor: A matrix containing the Taylor-like coefficients, where
% different orders are represented by different columns.
    [Nu, NC] = size(Cpoly);
    Ctaylor = zeros(Nu, NC);
    Ctaylor(:, 1) = Cpoly(:, 1);
    for polyi = 2:NC
        Ctaylor(:, 1:polyi) = Ctaylor(:, 1:polyi) + Cpoly(:, polyi)*Cp2t(polyi, 1:polyi);
    end
end
