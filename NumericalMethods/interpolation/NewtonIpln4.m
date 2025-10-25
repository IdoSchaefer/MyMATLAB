function result = NewtonIpln4(sp, fsp, resultp)
% The program performs a newton interpolation. It is intended also for
% vector functions. The program uses a domain of size 4 for numerical
% stability.
% sp: a vector of the sampling points.
% fsp: the vector of the values of the function in the sampling points.
% resultp: the point on which the desired function is to be evaluated.
    N = size(sp, 2);
    Dsize = max(sp) - min(sp);
    sp4 = 4/Dsize*sp;
    resultp4 = 4/Dsize*resultp;
    a = devdif(sp4, fsp);
    result = a(:, N);
    for spi = (N-1):-1:1
        result = a(:, spi) + (resultp4 - sp4(spi))*result;
    end
end