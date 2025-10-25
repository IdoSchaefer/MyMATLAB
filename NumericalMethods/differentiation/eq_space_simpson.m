function int = eq_space_simpson(f, a, b)
% The program computes the Simpson approximation of an integral, from
% equally spaced sampling points.
% The program assumes even N (odd number of sampling points, including the
% boundaries).
% Input:
% f: a vector of the function values in the equally spaced sampling points
% in the integration interval.
% a: the lower limit of integration
% b: the upper limit of integration
    N = length(f) - 1;
    % The space between neighbouring x points is h = (b - a)/N. 
    int = ((f(1) + 4*sum(f(2:2:N)) + 2*sum(f(3:2:N))+ f(N + 1)))*(b - a)/(3*N);
end