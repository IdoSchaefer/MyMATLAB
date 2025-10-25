function result = trig_fun(theta, A, B, s)
% The function returns the trigonometric function described in Fletcher and
% Powell, "A rapidly convergent descent method for minimization?", 1963.
% theta: A column vector which represents the variables (the alpha_j's in the paper)
% A, B: Matrices which represent the corresponding parameters in the paper
% s: A column vector which represents the column solution vector in the
% paper (the E_i's)
    result = sum((s - A*sin(theta) - B*cos(theta)).^2);
end