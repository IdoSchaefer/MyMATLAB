function grid = map_grid(eqfun, Ngrid)
% The function returns the grid points of a coordinate mapped into an
% equally spaced grid. It uses an analytical form of the mapping function
% from the original grid to the equally spaced grid. It is assumed that the
% potential has a center of inversion.
% eqfun: a function handle; the analytical mapping function from the 
% original variable of the grid to the equally spaced variable.
% Ngrid: the desired number of grid points. Assumed to be even.
% grid: the output grid in the original variable. The center
% of inversion is placed in 0, and its index is Ngrid/2+1.
    grid = zeros(Ngrid, 1);
    grid0i = floor(Ngrid/2) + 1;
    defaultop = optimset('fsolve');
    options = optimset(defaultop, 'Display', 'off');
    grid(grid0i) = fsolve(eqfun, 0, options);
    for gridi = 1:(Ngrid/2 - 1)
        grid(grid0i + gridi) = fsolve(@(x) eqfun(x) - gridi, grid(grid0i + gridi - 1), options);
    end
    grid((grid0i - 1):-1:2) = -grid((grid0i + 1):Ngrid);
    grid(1) = fsolve(@(x) eqfun(x) + Ngrid/2, grid(2), options);
end