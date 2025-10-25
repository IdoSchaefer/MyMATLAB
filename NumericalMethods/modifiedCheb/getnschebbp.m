function x = getnschebbp(Nx, alpha, beta, leftb, rightb)
% The function returns a nonsymmetric modified Chebyshev grid that includes the boundaries.
% Nx: The number of grid points minus 1.
% alpha: The streching parameter for the right boundary
% beta: The streching parameter for the left boundary
% leftb, rightb: The boundaries of the grid domain, [leftb rightb]
    xcheb = cos((0:Nx).'*pi/Nx);
    gtild1 = gtild(1, alpha, beta);
    gtildminus1 = gtild(-1, alpha, beta);
    a = 0.5*(gtild1 - gtildminus1);
    b = 0.5*(gtild1 + gtildminus1);
    x = ((gtild(xcheb, alpha, beta) - b)/a*(rightb - leftb) + leftb + rightb)/2;
end