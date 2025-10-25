function result = softrectfun1(x, a, b, stiffness_a, stiffness_b)
% The function cumputes a smooth rectangle function. May be non-symmetric,
% with different stiffness at the two boundaries.
    result = 0.5*(tanh(stiffness_a*(x - a)) - tanh(stiffness_b*(x - b)));
end