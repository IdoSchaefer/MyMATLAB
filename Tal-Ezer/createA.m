function A = createA(N)
    A = N^2*(diag(-2*ones(1, N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1));
end