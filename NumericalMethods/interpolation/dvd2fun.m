function result = dvd2fun(sp, polcoef, resultp)
% The program computes the Newton interpolation from the divided differences.
% It applies also for functions which return a vector. In such a case, the
% divided differences are represented as column vectors.
% sp: a vector of the sampling points. The last sampling point is
% unrequired, but can be nevertheless included in sp.
% polcoef: A matrix containing the vector divided differences in separate
% columns.
% resultp: A row vector; the points at which the desired function is to be evaluated.
    % The number of sampling points:
    N = size(polcoef, 2);
    % The degree of the interpolation polynomial is N-1.
    Nrp = length(resultp);
    result = repmat(polcoef(:, N), 1, Nrp);
    for spi = (N-1):-1:1
        result = polcoef(:, spi) + result.*(resultp - sp(spi));
    end
end