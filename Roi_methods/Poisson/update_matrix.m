function M = update_matrix(N, omega)
% The function creates an (N-1)^2 by (N-1)^2 matrix, which determines the
% update rule in an iterative smoothing procedure for solving the for the 
% Poisson equation in 2 dimensions. The boundary conditions are Phi = 0 on
% the boundaries.
% N: N+1 is the number of grid points, including the boundaries.
% omega: the smoothing parameter.
    I = eye(N-1);
    % The update matrix in one dimension:
    Mx = spdiags(ones(N - 1, 2), [-1 1], N - 1, N - 1);
    % The full update matrix is given by a direct product of the update
    % matrix in each dimension by the identity in the other dimension, 
    % with an addition of a diagonal term.
    M = (kron(Mx, I) + kron(I, Mx) + spdiags(2*omega*ones((N - 1)^2, 1), 0, (N - 1)^2, (N - 1)^2))/(2*(2 + omega));
end