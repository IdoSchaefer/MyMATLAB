function result = bFdirect(z, M)
    result = exp(z);
    for polyi = 1:M
        result = polyi*(result - 1)./z;
    end
end