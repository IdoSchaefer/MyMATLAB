function condition = terminate_reldif(dif_f, dif_x, fmin, gradmin, xmin, tolf, tolx)
% The function defines a termination condition for the quasiNewton program.
% The output is a boolean.
% dif_f: The difference of the function value from the previous iteration
% dix_x: A vector; the difference of the vector solution x from the previous iteration.
% fmin: The function value in the current iteration
% gradmin: The gradient in the current iteration
% xmin: The solution vector in the current iteration
% tolf: A tolerance parameter for the relative magnitude of the improvement
% in the function value
% tolx: A tolerance parameter for the relative magnitude of the improvement
% in the solution parameter
    condition = (abs(dif_f)/abs(fmin)<=tolf && norm(dif_x)/norm(xmin)<=tolx) || min(gradmin==0);
    % If min(gradmin==0) == true, it indicates that all the enteries of
    % gradmin are 0.
end