function result = softrectfun(x, a, b, stiffness)
% The function cumputes a smooth rectangle function.
    result = 0.5*(tanh(stiffness*(x - a)) - tanh(stiffness*(x - b)));
end