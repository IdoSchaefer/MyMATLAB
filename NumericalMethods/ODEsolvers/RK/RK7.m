function U = RK7(fderiv, t0tf, u0, dt, varargin)
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
    % The vector of the c coefficients, excluding the first term:
    %c = [1, 1, 2, 4, 6, 8, 10, 12]/12;
    c = [1/6, 1/3, 1/2, 2/11, 2/3, 6/7, 0, 1];
    % The K terms, which represent different approximations of du:
    K = zeros(dim, 9);
    % The vector of the b coefficients (weight of K terms):
%    weights = [41; 0; 0; 216; 27; 272; 27; 216; 41]/840;
    b8 = 77/1440;
    weights = [77/1440 - b8; 0; 0; 32/105; 1771561/6289920; 243/2560; 16807/74880; b8; 11/270];
    % The transposed alpha matrix, excluding the first row and the last column, which are 0:
%     alphaT = [1/12       0       0       0           0       0       0       0;
%               -5/6       11/12   0       0           0       0       0       0;
%               0          0       1/6     0           0       0       0       0;
%               157/9      -318/9  4/9     160/9       0       0       0       0;
%               -322/30    199/30  108/30  -131/30     0       0       0       0;
%               3158/45    -638/6  -69/6   314/6       0       157/45  0       0;
%               -53/14     38/7  -15/70  0           -65/72  0       29/90   0;
%               56/25      849/42  -833/42 -156/42     -39/45  149/32  -125/45 27/25].';
    alphaT = [1/6                       0       0           0               0                   0               0                   0;
              0                         1/3     0           0               0                   0               0                   0;
              1/8                       0       3/8         0               0                   0               0                   0;
              148/1331                  0       150/1331    -56/1331        0                   0               0                   0;
              -404/243                  0       -170/27     4024/1701       10648/1701          0               0                   0;
              2466/2401                 0       1242/343    -19176/16807    -51909/16807        1053/2401       0                   0;
              1/(576*b8)                0       0           1/(105*b8)      -1331/(279552*b8)   -9/(1024*b8)    343/(149760*b8)     0;
              -71/32 - 270*b8/11        0       -195/22     32/7            29403/3584          -729/512        1029/1408           270*b8/11].';
    for ti = 1:Nt
        K(:, 1) = fderiv(t, U(:, ti), varargin{:})*dt;
        for Ki = 2:9
            K(:, Ki) = fderiv(t + c(Ki - 1)*dt, U(:, ti) + K(:, 1:(Ki - 1))*alphaT(1:(Ki - 1), Ki - 1), varargin{:})*dt;
        end
        U(:, ti + 1) = U(:, ti) + K*weights;
        t = t + dt;
    end   
end