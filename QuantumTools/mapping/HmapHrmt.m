function H = HmapHrmt(x, p, Vvec, Jvec, DJvec, D2Jvec, m)
% The function computes the "Hamiltonian" for a mapped grid, where the
% eigenvectors are J^(1/2)*psi, where psi is the quantum vector state, and 
% J is the Jacobian. In this form, the "Hamiltonian" is Hermitian.
% x: the x grid in the equally spaced grid.
% p: the p grid in the equally spaced grid.
% Vvec: a vector that contains the potential energy values in all grid
% points. It should be computed in the original grid r.
% Jvec: a vector that contains the Jacobian values in all grid points. It 
% should be computed in the original grid r.
    if nargin < 7
        m = 1;
    end
    Nx = length(x);
    invJsq = spdiags(1./(Jvec.^2), 0, Nx, Nx);
    Psqx = Nx*ifft(ifft(diag(p.^2))')';
    H = diag(Vvec + (7*DJvec.^2./(2*Jvec) - D2Jvec)./(4*m*Jvec.^3)) + (Psqx*invJsq + invJsq*Psqx)/(4*m);
end