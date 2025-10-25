function exval_mu = evmu(psi, mu)
% The function computes the time-dependent expectation value of the dipole
% moment exval_mu from the time-dependent state vector psi and the dipole
% moment mu.
% psi: a matrix; the state vector at all time points, where different
% time-points are represented by separate columns. Represented in the x
% domain.
% mu: a column vector; the x dependent dipole moment function, represented
% in the x domain
    exval_mu = sum(psi.*conj(psi).*mu);
end