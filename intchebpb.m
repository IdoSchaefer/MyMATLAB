function result = intchebpb(fv, lengthD)
% The function computes the integral on a function sampled in Chebyshev
% points that include the boundary.
% If you have to compute more than one function sampled in the same number
% of Chebyshev points, use the program: chebweights, to get the weights of
% integration (which are independent of the function).
% fv: the values of the function in the Chebyshev points.
% lengthD: the length of the interval of the integration. 
    dim = size(fv);
    if dim(1) == 1
        N = dim(2);
        n_even = 0:2:(N-1);
    else
        N = dim(1);
        n_even = (0:2:(N-1)).';
    end
% Getting the Chebyshev coefficients:
    c = chebcbv(fv);
% The integral on the Chebyshev polynomials can be performed analytically,
% and is found to be:
% 0             for odd n
% -2/(n^2 - 1)  for even n
% The integral on the domain [-1 1] is multiplied by a factor of: lengthD/2,
% to get the result in the correct domain of integration.
    result = -sum(c(n_even + 1)./(n_even.^2 - 1))*lengthD;
end