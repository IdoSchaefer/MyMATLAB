function a = newtoncoef(x, f)
    N = length(x);
    a = zeros(1, N);
    last_diag = zeros(1, N);
    new_diag = zeros(1, N);
    for i = 1:N
        new_diag(1) = f(i);
        for j = 2:i
            new_diag(j) = (new_diag(j-1) - last_diag(j-1))/(x(i) - x(i-j+1));
        end
        a(i) = new_diag(i);
        last_diag = new_diag;
    end
end
