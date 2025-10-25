function mE = evHt(psi, K, V, field, mu)
% The function computes the Hamiltonian expectation value in each time step
% of psi, including a time dependent field.
% V: the potential energy. It's a vector in the x domain.
% K: the kinetic energy vector in the p domain. Has to be computed in the
% following way:
%     p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
%     p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
%     K = p.^2/2;
% field: a vector of the time dependent field.
% mu: the vector of the dipole function.
    Nt = size(psi, 2);
    mE = zeros(1, Nt);
    if nargin<4
        field = zeros(1, Nt);
        mu = 0;
    end
    for ti = 1:Nt
        mE(ti) = psi(:, ti)'*Hpsi(K, V - field(ti)*mu, psi(:, ti));
    end
end