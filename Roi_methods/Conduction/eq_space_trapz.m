function int = eq_space_trapz(f, a, b)
% The program computes the trapeze approximation of an integral, from
% equally spaced sampling points.
% Input:
% f: a vector of the function values in the equally spaced sampling points
% in the integration interval.
% a: the lower limit of integration
% b: the upper limit of integration
    N = length(f) - 1;
    % The space between neighbouring x points is h = (b - a)/N. 
    int = ((f(1) + f(N + 1))/2 + sum(f(2:N)))*(b - a)/N;
end