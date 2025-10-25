function result = NewtonIpln4b(sp, fsp, resultp)
% The program performs a newton interpolation. It is intended also for
% vector functions. The interpolation can be performed in more than one point.
% The program uses a domain of size 4 for numerical
% stability.
% sp: a vector of the sampling points.
% fsp: the vector of the values of the function in the sampling points.
% resultp: a vector of the points on which the desired function is to be evaluated.
    Nsp = length(sp);
    Nresult = length(resultp);
    Dsize = max(sp) - min(sp);
    % The order of points is mixed in a random way, for numerical stability:
    [sp, spi] = rand_order(sp);
    fsp = fsp(:, spi);
    sp4 = 4/Dsize*sp;
    resultp4 = 4/Dsize*resultp;
    a = devdif(sp4, fsp);
    result = a(:, Nsp)*ones(1, Nresult);
    resultp4M = spdiags(resultp4.', 0, Nresult, Nresult);
    for spi = (Nsp-1):-1:1
        result = a(:, spi)*ones(1, Nresult) + result*resultp4M - sp4(spi)*result;
    end
end