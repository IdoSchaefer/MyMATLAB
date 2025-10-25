function [Vnew, Vnewx] = newguess0lb3(Vold, Nnew)
% The function returns a new potential guess for a denser grid, based on a
% result from a coarse grid. Intended for percosV0lb_grad3.
% Vold: The optimized potential in the coarse grid.
% Nnew: The number of points in the new grid minus 1.
    dim = size(Vold);
    if dim(1) == 1
        Nold = dim(2)/2;
        Vold = Vold.';
    else
        Nold = dim(1)/2;
    end
    Vnew = [Vold(1:Nold); -(Vold(1)/2 + sum(Vold(2:Nold))); zeros(Nnew - Nold - 1, 1);...
        Vold((Nold + 1):2*Nold); -(Vold(Nold + 1)/2 + sum(Vold((Nold + 2):2*Nold))); zeros(Nnew - Nold - 1, 1)]*sqrt(Nnew/Nold);
    Vnewx = dctI0lb(Vnew(1:Nnew)) - 1i*dctI0lb(Vnew((Nnew + 1):2*Nnew)).^2;
    if dim(1) == 1
        Vnew = Vnew.';
        Vnewx = Vnewx.';
    end
end