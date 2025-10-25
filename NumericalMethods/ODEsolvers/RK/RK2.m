function U = RK2(fderiv, t0tf, u0, dt, varargin)
% fderiv: a function handle of the form @(t, u, varargin). The value of the function
% is the derivative of u with respect to t.
% t0tf:  a vector that contains initial and final time: [t0 tf]
% u0:  a vector of the initial values of the result vector u.
% dt: The time step.
% U:  The result u vectors at all the times.
    t0 = t0tf(1);
    tf = t0tf(2);
    Nt = round(abs((tf - t0)/dt));
    dim = length(u0);
    U = zeros(dim, Nt + 1);
    U(:, 1) = u0;
    t = t0;
    for ti = 1:Nt
        U(:, ti + 1) = U(:, ti) + fderiv(t + dt/2, U(:, ti) + fderiv(t, U(:, ti))*dt/2, varargin{:})*dt;
        t = t + dt;
    end
end