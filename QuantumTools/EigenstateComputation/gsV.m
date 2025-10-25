function [fi0, E0, x, E, P, H] = gsV(Vf, xdomain, Nx, m)
% The function returns the ground state function and energy.
% Vf: The potential. A function handle of the form: @(x), where x is a
% vector. It is possible to insert the potential vector itself instead.
% xdomain: the domain of x. xdomain = [minx maxx].
% Nx: The number of x points.
% m: the mass
    if nargin < 4
        m = 1;
    end
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
%    V = zeros(Nx, 1);
    % Making a vector of the potential energy at all x:
    if length(Vf) == 1
        % If Vf is a function handle:
        V = Vf(x);
    else
        % If Vf is a vector:
        V = Vf;
    end
%    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
%    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    p = (-pi/dx:(2*pi/xdlength):(pi/dx - 2*pi/xdlength)).';
    % Making a vector of the kinetic energy at all p, for m=1:
    K = p.^2/(2*m);
    % The hamiltonian in the x domain:
    H = diag(V) + p2x(diag(K));
    [P, D] = eig(H);
    E = diag(D);
%    [E, orderE] = sort(E);
    [Ereal, orderE] = sort(real(E));
    E = E(orderE);
    P = P(:, orderE);
    E0 = E(1);
    fi0 = P(:, 1);
%     [E0, nmin] = min(E);
%     fi0 = P(:, nmin);
end