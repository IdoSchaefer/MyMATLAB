function u = RK7uf(fderiv, t0tf, u0, dt, varargin)
% fderiv: a function handle of the form @(t, u, varargin). The value of the function
% is the derivative of u with respect to t.
% t0tf:  a vector that contains initial and final time: [t0 tf]
% u0:  a vector of the initial values of the result vector u.
% dt: The time step.
% u:  The result u vector at the final time.
    t0 = t0tf(1);
    tf = t0tf(2);
    Nt = round(abs((tf - t0)/dt));
    dim = length(u0);
    u = u0;
    t = t0;
    % The vector of the c coefficients, excluding the first term:
    c = [1/6, 1/3, 1/2, 2/11, 2/3, 6/7, 0, 1];
    % The K terms, which represent different approximations of du:
    K = zeros(dim, 9);
    % The vector of the b coefficients (weight of K terms):
    b8 = 77/1440;
    weights = [77/1440 - b8; 0; 0; 32/105; 1771561/6289920; 243/2560; 16807/74880; b8; 11/270];
    % The transposed alpha matrix, excluding the first row and the last column, which are 0:
    alphaT = [1/6                       0       0           0               0                   0               0                   0;
              0                         1/3     0           0               0                   0               0                   0;
              1/8                       0       3/8         0               0                   0               0                   0;
              148/1331                  0       150/1331    -56/1331        0                   0               0                   0;
              -404/243                  0       -170/27     4024/1701       10648/1701          0               0                   0;
              2466/2401                 0       1242/343    -19176/16807    -51909/16807        1053/2401       0                   0;
              1/(576*b8)                0       0           1/(105*b8)      -1331/(279552*b8)   -9/(1024*b8)    343/(149760*b8)     0;
              -71/32 - 270*b8/11        0       -195/22     32/7            29403/3584          -729/512        1029/1408           270*b8/11].';
    for ti = 1:Nt
        K(:, 1) = fderiv(t, u, varargin{:})*dt;
        for Ki = 2:9
            K(:, Ki) = fderiv(t + c(Ki - 1)*dt, u + K(:, 1:(Ki - 1))*alphaT(1:(Ki - 1), Ki - 1), varargin{:})*dt;
        end
        u = u + K*weights;
        t = t + dt;
    end   
end