function L = Laplacian_matrix(N, h)
% The function creates a finite difference approximation to the Laplacian 
% operator in 2 dimensions. The boundary conditions are Phi = 0 on
% the boundaries.
% N+1 is the number of grid points, including the boundaries.
% h: The spacing between grid points
    I = eye(N-1);
    % The Laplacian matrix in one dimension:
    ones_vec = ones(N - 1, 1);
    Lx = spdiags([ones_vec, -2*ones_vec, ones_vec], -1:1, N - 1, N - 1)/h^2;
    % The full Laplacian matrix is given by a direct product of the Laplacian
    % matrix in each dimension by the identity in the other dimension: 
    L = kron(Lx, I) + kron(I, Lx);
end