function [Vnew, Vnewx] = newguess0lb_sp(Vold, Nnew)
% The function returns a new potential guess for a denser grid, based on a
% result from a coarse grid. Intended for percosV0lb_grad. A spline
% interpolation is used.
% Vold: The optimized potential in the coarse grid.
% Nnew: The number of points in the new grid minus 1.
    dim = size(Vold);
    if dim(1) == 1
        Nold = dim(2)/2;
        Vold = Vold.';
    else
        Nold = dim(1)/2;
    end
    Voldx = dctI0lb(Vold(1:Nold)) + 1i*dctI0lb(Vold((Nold + 1):2*Nold));
    Vnewx = spline((0:Nold).', Voldx, (0:Nold/Nnew:Nold).');
    Vnew = Vx2Vk0lb([real(Vnewx(2:(Nnew + 1))); imag(Vnewx(2:(Nnew + 1)))]);
    if dim(1) == 1
        Vnew = Vnew.';
        Vnewx = Vnewx.';
    end
end