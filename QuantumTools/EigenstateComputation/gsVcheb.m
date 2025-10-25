function [fi0, E0, x, E, P, H] = gsVcheb(Vf, xdomain, Nx, m)
% The function returns the ground state function and energy.
% Vf: The potential. A function handle of the form: @(x), where x is a
% vector.
% xdomain: the domain of x. xdomain = [minx maxx].
% Nx: The number of x points.
% m: the mass
    if nargin < 4
        m = 1;
    end
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    theta = ((1:Nx).'*2 - 1)*pi/(2*Nx);
    xcheb = -cos(theta);
    x = 0.5*(xcheb*xdlength + min_x + max_x);
    % Making a vector of the potential energy at all x:
    V = Vf(x);
    % The hamiltonian in the x domain:
    D2 = D2chebMat(Nx, xdlength);
    H = diag(V) - D2/(2*m);
    %metric = spdiags(sin(theta.')*pi/Nx*xdlength/2, 0, Nx, Nx);
    %metric = diag(sin(theta.')*pi/Nx*xdlength/2);
    %[P, D] = eig(metric*H, metric);
    [P, D] = eig(H);
    E = diag(D);
    [E, orderE] = sort(E);
    P = P(:, orderE);
    E0 = E(1);
    fi0 = P(:, 1);
end