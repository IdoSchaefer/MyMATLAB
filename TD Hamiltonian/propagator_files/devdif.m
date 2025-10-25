function polcoef = devdif(x, fx)
    [dim, Npoints] = size(fx);
    polcoef = zeros(dim, Npoints);
    diagonal = zeros(dim, Npoints);
    polcoef(:, 1) = fx(:, 1);
    diagonal(:, 1) = fx(:, 1);
    for coefi = 2:Npoints
        diagonal(:, coefi) = fx(:, coefi);
        for dtermi = coefi-1:-1:1
            diagonal(:, dtermi) = (diagonal(:, dtermi + 1) - diagonal(:, dtermi))./(x(coefi) - x(dtermi));
% The value of diagonal(:, dtermi) is from the previous diagonal, above diagonal(:, dterm + 1).
        end
        polcoef(:, coefi) = diagonal(:, 1);
    end
end