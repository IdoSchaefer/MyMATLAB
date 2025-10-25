function x = getmodchebbp(Nx, p, leftb, rightb)
% The function returns a modified Chebyshev grid that includes the boundaries.
% Nx: The number of grid points minus 1.
% p: The streching parameter from the new article, 0<p<pi/2
% leftb, rightb: The boundaries of the grid domain, [leftb rightb]
    xcheb = cos((0:Nx).'*pi/Nx);
    x = (1/p*asin(xcheb*sin(p))*(rightb - leftb) + leftb + rightb)/2;
end