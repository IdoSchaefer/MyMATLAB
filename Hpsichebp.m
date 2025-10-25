function result = Hpsichebp(V, psi, domainL, m)
% The result is the operation of the Hamiltonian on the wave function
% psi, sampled in the Chebyshev points.
% V: The potential energy. It's a vector in the x domain.
% domainL: The length of the x domain.
% m: The mass
    result = V.*psi - D2cheb(psi, domainL)/(2*m);
end