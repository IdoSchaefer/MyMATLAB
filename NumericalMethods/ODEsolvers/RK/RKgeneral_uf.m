function u = RKgeneral_uf(fderiv, t0tf, u0, dt, c, weights, alphaT, varargin)
% The function computes the RK solution of and ODE system; applies for the
% general RK scheme with values of the Butcher table as an input.
% Input:
% fderiv: a function handle of the form @(t, u, varargin). The value of the function
% is the derivative of u with respect to t.
% t0tf:  a vector that contains initial and final time: [t0 tf]
% u0:  a vector of the initial values of the result vector u.
% dt: The time step
% c: The vector of the c coefficients, excluding the first term, which is 0.
% weight: Column vector; the vector of the b coefficients (weight of K terms)
% alphaT: The transposed alpha matrix, excluding the first row and the last
% column of alpha, which are 0.
% Output:
% U:  The result u vectors at all the times.
    t0 = t0tf(1);
    tf = t0tf(2);
    Nt = round(abs((tf - t0)/dt));
    dim = length(u0);
    u = u0;
    t = t0;
    NK = length(weights);
    % The K terms, which represent different approximations of du:
    K = zeros(dim, NK);
    for ti = 1:Nt
        K(:, 1) = fderiv(t, u, varargin{:})*dt;
        for Ki = 2:9
            K(:, Ki) = fderiv(t + c(Ki - 1)*dt, u + K(:, 1:(Ki - 1))*alphaT(1:(Ki - 1), Ki - 1), varargin{:})*dt;
        end
        u = u + K*weights;
        t = t + dt;
    end   
end