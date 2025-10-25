function [r, Nr, Q] = map_grid1(eqfun, jfun, maxr)
%function [r, Nr, Q] = map_grid1(eqfun, maxr)
% The function returns the grid points of a coordinate mapped into an
% equally spaced grid. It uses an analytical form of the mapping function
% from the original grid to the equally spaced grid. It is assumed that the
% potential has a center of inversion.
% Input:
% eqfun: a function handle; the analytical mapping function from the 
% original variable of the grid to the equally spaced variable.
% jfun: a function handle; the jacobian of the mapping function.
% maxr: the maximal allowed value of the original r grid.
% Output:
% r: the output grid in the original variable. The center
% of inversion is placed in 0 in the equally spaced grid, and its index is Nr/2+1.
% Nr: the number of grid points.
% Q: the equidistant grid.
    defaultop = optimset('fsolve');
    options = optimset(defaultop, 'Display', 'off', 'jacobian', 'on');
%    options = optimset(defaultop, 'Display', 'off');
%     r0 = fsolve(eqfun, 0, options);
%     r1 = fsolve(@(x) eqfun(x) - 1, r0, options);
    r0 = fsolve(@(x) solver_input(x, 0), 0, options);
    r1 = fsolve(@(x) solver_input(x, 1), r0, options);
    dr = r1 - r0;
    r = zeros(ceil(maxr/dr), 1);
    r(1) = r0;
    ri = 2;
    r(ri) = r1;
%    dr = r(2) - r(1);
    while r(ri)<=maxr
        ri = ri + 1;
%        drprevious = dr;
        dr = r(ri - 1) - r(ri - 2);
%         d2r = dr - drprevious;
%        r(ri) = fsolve(@(x) eqfun(x) - ri + 1, r(ri - 1) + dr, options);
        r(ri) = fsolve(@(x) solver_input(x, ri - 1), r(ri - 1) + dr, options);
    end
    r = [-r((ri - 1):-1:2); r(1:(ri-2))];
    Nr = length(r);
    Q = (-Nr/2):(Nr/2 - 1);

%%%% Nested function %%%%
    function [funval, jval] = solver_input(r, Q)
        funval = eqfun(r) - Q;
        jval = jfun(r);
    end    
end