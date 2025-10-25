function result = Pchebpb(x, N)
    result = (x.^2 - 1);
    if N > 2
        U0 = ones(size(x));
        U1 = 2*x;
    end
    for Ui = 2:(N - 2)
        U2 = 2*x.*U1 - U0;
        U0 = U1;
        U1 = U2;
    end
    if N>=4
        result = result.*U2/2^(N - 2);
    elseif N == 3
        result = result.*U1/2^(N - 2);
    end
end
