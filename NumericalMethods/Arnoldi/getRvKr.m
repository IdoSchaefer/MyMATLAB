function Rv = getRvKr(Hessenberg, v, samplingp, Nkr, capacity)
% For the Newton approximation of a function of matrix which multiplies a
% vector:
% f(A)v \approx \sum_{n=0}^Nkr a_n*R_n(A)v, 
% the function computes the R_n(A)v vectors represented in the Krylov
% space of dimension Nkr+1, where the R_n(z)'s are the Newton basis 
% polynomials, with samplingp as the sampling points. The Newton
% approximation is performed in a space of capacity 1.
% Input:
% Hessenberg: The extended Hessenberg matrix of the problem
% v: Column vector; defined mathematically above.
% samplingp: The sampling points; should be the eigenvalues of the
% Hessenberg matrix.
% Nkr: The dimension of the Krylov space which is actually used for the Arnoldi approxiamtion.
% capacity: The capacity of the approximation domain (see get_capacity.m).
    Rv = zeros(Nkr + 1, Nkr + 1);
    % Rv(:, 1) is v in the Krylov space.
    Rv(1, 1) = norm(v);
    for spi = 1:Nkr
        % Rv(:, spi) belongs to a Krylov space of dimension spi, and
        % the terms with row indices larger than spi vanish.
        Rv(1:(spi + 1), spi + 1) = (Hessenberg(1:(spi + 1), 1:spi)*Rv(1:spi, spi) - samplingp(spi)*Rv(1:(spi + 1), spi))/capacity;
    end
end
