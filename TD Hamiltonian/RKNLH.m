function U = RKNLH(H0, Vt, u0, x, T, dt)
% The program solves The Schrodinger equation, with a non-linear Hamiltonian of the
% form: H = H0 + Vt. Vt is the disturbance. It's an analytic function of
% the variable of u0 (x in the common case). It depentds on time explicitly
% and implicitly, through the dependence on the state vector u.
% H0 is a matrix. Vt is a function handle of the form: @(u, x, t). u0 is the
% initial state vector. x is the grid of u0. T is the total time. dt is the
% time difference between the time points that will be returned.
% The columns of the solution, U, represent different time points.
    if nargin<6
        dt = 0.1;
    end
    szx = size(x);
    if szx(2)>1
        x = x.';
    end
    sol = ode45(@NLH, [0 T], u0, [], H0, Vt, x);
    U = deval(sol, 0:dt:T);
end