function [Vnew, Vnewx] = newguess0lb(Vold, Nnew)
% The function returns a new potential guess for a denser grid, based on a
% result from a coarse grid. Intended for percosV0lb_grad.
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
    ReVnewx = dctI0lb(Vnew(1:Nnew));
    ImVnewx = dctI0lb(Vnew((Nnew + 1):2*Nnew));
    if max(ImVnewx)<=0
        Vnewx = ReVnewx + 1i*ImVnewx;
    else
    % If there are positive imaginary values, a spline interpolation is
    % used instead of the cosine interpolation:
        display('The cosine interpolation results in positive imaginary values. Using a spline interpolation instead for the imaginary part.')
        ImVoldx = dctI0lb(Vold((Nold + 1):2*Nold));
        ImVnewx = spline((0:Nold).', ImVoldx, (0:Nold/Nnew:Nold).');
        Vnewx = ReVnewx + 1i*ImVnewx;
        Vnew = Vx2Vk0lb([real(Vnewx(2:(Nnew + 1))); imag(Vnewx(2:(Nnew + 1)))]);
    end
    if dim(1) == 1
        Vnew = Vnew.';
        Vnewx = Vnewx.';
    end
end