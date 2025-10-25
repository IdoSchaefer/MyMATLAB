function [fi0, E0, x, E, P, H] = gsVgen(Hfun, xdomain, Nx, varargin)
% The function returns the ground state function and energy, for a
% generalized form of the "Hamiltonian", which includes an arbitrary
% dependence on x and p.
% Hfun: a function that creates the Hamiltonian matrix from the x and p
% grids.
% xdomain: the domain of x. xdomain = [minx maxx].
% Nx: The number of x points.
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    % The hamiltonian in the x domain:
    H = Hfun(x, p, varargin{:});
    [P, D] = eig(H);
    E = diag(D);
    [E, orderE] = sort(E);
    P = P(:, orderE);
    E0 = E(1);
    fi0 = P(:, 1);
end