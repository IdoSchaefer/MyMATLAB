function plot_steps(x, y)
    szx = size(x);
    szy = size(y);
    if szx(1) == 1
        Nsteps = szx(2) - 1;
        extended_x = [x(1), kron(x(2:Nsteps), [1, 1]), x(Nsteps + 1)];
    else
        Nsteps = szx(1) - 1;
        extended_x = [x(1); kron(x(2:Nsteps), [1; 1]); x(Nsteps + 1)];
    end
    if szy(1) == Nsteps
        extended_y = kron(y, [1; 1]);
    elseif szy(2) == Nsteps
        extended_y = kron(y, [1, 1]);
    else
        fprintf('\nError: Input dimension mismatch.\n')
        return
    end
    plot(extended_x, extended_y)
end