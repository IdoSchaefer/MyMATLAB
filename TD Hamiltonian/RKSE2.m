function U = RKSE2(Vt, u0, xdomain, T, dt)
% The program solves The Schrodinger equation, with an Hamiltonian of the
% form: H = T + Vt. Vt is the potential energy. It's an analytic function of
% the variable of u0 (x in the common case).
% Vt is a function handle of the form: @(x, t). u0 is the
% initial state vector. xdomain is the domain of x. xdomain = [minx maxx].
% T is the total time. dt is the
% time difference between the time points that will be returned.
% The columns of the solution, U, represent different time points.
    Nx = length(u0);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
    if nargin<5
        dt = 0.1;
    end
    sol = ode45(@SE2, [0 T], u0, [], K, Vt, x);
    U = deval(sol, 0:dt:T);
end