function M = op2mat(op, N)
% The function finds the matrix representation M of a linear operator op.
% op is a function handle.
% N: the dimension of the matrix.
    en = zeros(N, 1);
    M = zeros(N, N);
    for ni = 1:N
        en(ni) = 1;
        M(:, ni) = op(en);
        en(ni) = 0;
    end
end