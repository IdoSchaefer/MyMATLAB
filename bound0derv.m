function result = bound0derv(v, xdomain, a, b, stiffness, Nxfactor)
% The function computes a variant of the signal v with gradually
% decaying first derivative at the edges.
% This program can be made more accurate by using interpolation
% Input:
% v: A vector.
% x: The x grid
% a: The midpoint of the region of the decay of the derivative near the
% left boundary.
% b: The midpoint of the region of the decay of the derivative near the
% right boundary.
% stiffness: The stiffness of the smoothing function.
% Output:
% result: A vector which contains the resulting function values in the x
% points.
    Nx = length(v) - 1;
    dx = (xdomain(2) - xdomain(1))/Nx;
    x = xdomain(1):dx:xdomain(2);
    intdx = dx/Nxfactor;
    intx = xdomain(1):intdx:xdomain(2);
    if size(v, 2) == 1
        x = x.';
        intx = intx.';
    end
    intv = spline(x, v, intx);
    integrand = intv.*0.5*stiffness.*(-sech(stiffness*(intx - a)).^2 + sech(stiffness*(intx - b)).^2);
    result = v.*softrectfun(x, a, b, stiffness);
    integx = 0;
    for xi = 2:(Nx + 1)
        intxi = ((xi - 2)*Nxfactor + 1):((xi - 1)*Nxfactor + 1);
        integx = integx + (0.5*(integrand(intxi(1)) + integrand(intxi(Nxfactor + 1))) + sum(integrand(intxi(2:(Nxfactor)))))*intdx;
        result(xi) = result(xi) + integx;
    end
    midpointi = round(Nx/2);
    result = result - result(midpointi) + v(midpointi);
end