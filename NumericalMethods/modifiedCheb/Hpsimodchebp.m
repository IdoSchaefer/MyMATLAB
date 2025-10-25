function result = Hpsimodchebp(V, psi, domainL, D2coef, Dcoef, m)
% The result is the operation of the Hamiltonian on the wave function
% psi, sampled in the modified Chebyshev points.
% V: The potential energy. It's a vector in the x domain.
% domainL: The length of the x domain.
% m: The mass
    result = V.*psi - (D2coef.*D2chebb(psi, domainL) + Dcoef.*Dchebb(psi, domainL))/(2*m);
end