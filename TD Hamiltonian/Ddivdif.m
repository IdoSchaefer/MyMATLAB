function result = Ddivdif(x, Dfx_N, diagonal)
    [~, Npoints] = size(diagonal);
    result = Dfx_N;
    for xi = (Npoints - 1):-1:1
        result = (result - diagonal(xi))/(x(Npoints) - x(xi));
    end
end