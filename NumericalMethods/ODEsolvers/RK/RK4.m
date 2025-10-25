function U = RK4(fderiv, t0tf, u0, dt, varargin)
% fderiv:  a function handle. The value of the function is the derivative of u with respect to t.
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
        ut = U(:, ti);
 % The k's are estimations of du. The result is composed of a weighted
 % average of the 4 k's.
        k1 = feval(fderiv, t, ut, varargin{:})*dt;
        k2 = feval(fderiv, t + 0.5*dt, ut + 0.5*k1, varargin{:})*dt;
        k3 = feval(fderiv, t + 0.5*dt, ut + 0.5*k2, varargin{:})*dt;
        k4 = feval(fderiv, t + dt, ut + k3, varargin{:})*dt;
        U(:, ti + 1) = U(:, ti) + (k1 + 2*k2 + 2*k3 + k4)/6;
        t = t + dt;
    end
end