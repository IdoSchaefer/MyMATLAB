function Dv = DchebbMap(v, r)
% The function returns the derivative of a periodic function sampled at the 
% points r, using mapped Chebyshev differentiation.
% v: A vector that contains the function values. It is assumed to satisfy the
% cosine series boundary conditions (zero derivative).
% r: A vector that contains the sampling points.
    Dv = Dchebb(v, 2)./Dchebb(r, 2);
end