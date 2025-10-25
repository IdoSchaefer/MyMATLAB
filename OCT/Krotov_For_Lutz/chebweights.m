function w = chebweights(N, lengthD)
% The program computes the weights of the function values sampled on
% Chebyshev points that include the boundary, used for integration.
% N:    The number of Chebyshev points.
% lengthD:  The length of the domain.
% Exact for a polynomial of degree N-1.
    integT = zeros(1, N);
    n_even = (0:2:(N-1));
    integT(n_even + 1) = -1./(n_even.^2 - 1);
    w = lengthD*chebcbv(integT);
end