function result = bound0der(fun, x, a, b, stiffness)
% The function computes a variant of the function fun with gradually
% decaying first derivative at the edges. This function computes only an
% approximation. The exact calculation is performed in bound0der1.
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
    fx = fun(x);
    result = fx + softstepfun(a - x, stiffness).*(fun(a) - fx)...
                  + softstepfun(x - b, stiffness).*(fun(b) - fx);
end