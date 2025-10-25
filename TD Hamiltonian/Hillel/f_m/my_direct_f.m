function result = my_direct_f(z, m)
    result = exp(z);
    for polyi = 1:m
        result = polyi*(result - 1)./z;
    end
end