function result = bound0der1(fun, x, a, b, stiffness)
% The function computes a variant of the function fun with gradually
% decaying first derivative at the edges.
% Input:
% fun: A function handle.
% x: The x grid
% a: The midpoint of the region of the decay of the derivative near the
% left boundary.
% b: The midpoint of the region of the decay of the derivative near the
% right boundary.
% stiffness: The stiffness of the smoothing function.
% Output:
% result: A vector which contains the resulting function values in the x
% points.
    fintegrand = @(x) fun(x)*0.5*stiffness.*(-sech(stiffness*(x - a)).^2 + sech(stiffness*(x - b)).^2);
    Nx = length(x);
    result = fun(x).*softrectfun(x, a, b, stiffness);
    integx = 0;
    for xi = 2:Nx
        integx = integx + integral(fintegrand, x(xi - 1), x(xi));
        %integx = integx + quad(fintegrand, x(xi - 1), x(xi));
        result(xi) = result(xi) + integx;
    end
    midpointi = round(Nx/2);
    result = result - result(midpointi) + fun(x(midpointi));
end