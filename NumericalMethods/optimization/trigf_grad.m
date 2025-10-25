function result = trigf_grad(theta, A, B, s)
% The function returns the gradient of the trigonometric function described in Fletcher and
% Powell, "A rapidly convergent descent method for minimization?", 1963.
% theta: A column vector which represents the variables (the alpha_j's in the paper)
% A, B: Matrices which represent the corresponding parameters in the paper
% s: A column vector which represents the column solution vector in the
% paper (the E_i's)
    v = s - A*sin(theta) - B*cos(theta);
    diagv = diag(v);
    result = 2*((sum(diagv*B).').*sin(theta) - (sum(diagv*A).').*cos(theta));
end