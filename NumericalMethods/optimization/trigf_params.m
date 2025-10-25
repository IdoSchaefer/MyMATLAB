function [A, B, s, theta_min] = trigf_params(N)
% The function generates the parameters for the function trig_fun,
% according to the paper of Fletcher and Powell.
    A = (rand(N) - 0.5)*200;
    B = (rand(N) - 0.5)*200;
    theta_min = (rand(N, 1) - 0.5)*2*pi;
    s = A*sin(theta_min) + B*cos(theta_min);
end