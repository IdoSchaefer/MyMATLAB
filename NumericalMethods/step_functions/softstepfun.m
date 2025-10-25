function result = softstepfun(x, stiffness)
% The function cumputes a smooth step function.
    result = 0.5*(tanh(stiffness*x) + 1);
end