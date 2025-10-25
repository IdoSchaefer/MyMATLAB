function result = NewtonIpln(sp, fsp, resultp)
% The program performs a newton interpolation. It is intended also for
% vector functions.
% sp: a vector of the sampling points.
% fsp: the vector of the values of the function in the sampling points.
% resultp: the point on which the desired function is to be evaluated.
    N = size(sp, 2);
    a = devdif(sp, fsp);
    result = a(:, N);
    for spi = (N-1):-1:1
        result = a(:, spi) + (resultp - sp(spi))*result;
    end
end