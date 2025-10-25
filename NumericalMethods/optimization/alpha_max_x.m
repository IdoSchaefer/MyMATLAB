function max_alpha = alpha_max_x(x0, direction, max_x)
% The function returns the maximal allowed alpha value for the line-search
% procedure line_search1, for a given maximal allowed value of the 
% magnitude of the optimization variables. It is possible to set a uniform
% value for all the variables (the infinity norm of x), or to specify the
% value for each variable.
% x0: The initial point of the search, i.e. the solution from the previous
% iteration
% direction: The direction of search
% max_x: The value of the maximal allowed magnitude of x for each of the x
% components; if all the variables are restricted to a uniform maximal
% magnitude, max_x may be a positive scalar. If there are different maximal
% magnitudes for the different variables, max_x is a vector of the
% dimension of x0.
    max_alpha = min((max_x - sign(direction).*x0)./abs(direction));
end