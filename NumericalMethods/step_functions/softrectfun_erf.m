function result = softrectfun_erf(x, a, b, sigma)
% The function cumputes a smooth rectangle function.
    result = 0.5*(erf((x - a)/(sqrt(2)*sigma)) - erf((x - b)/(sqrt(2)*sigma)));
end