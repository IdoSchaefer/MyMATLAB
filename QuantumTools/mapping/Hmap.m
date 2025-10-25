function H = Hmap(x, p, Vvec, Jvec, m)
% The function computes the Hamiltonian for a mapped grid.
% x: the x grid in the equally spaced grid.
% p: the p grid in the equally spaced grid.
% Vvec: a vector that contains the potential energy values in all grid
% points. It should be computed in the original grid r.
% Jvec: a vector that contains the Jacobian values in all grid points. It 
% should be computed in the original grid r.
    if nargin < 5
        m = 1;
    end
    Nx = length(x);
    invJ = spdiags(1./Jvec, 0, Nx, Nx);
    Px = Nx*ifft(ifft(diag(p))')';
    Pr = invJ*Px;
    H = diag(Vvec) + Pr^2/(2*m);
end