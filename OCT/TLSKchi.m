function chi = TLSKchi(E, chiT, psi, Epenal, T, Nt)
% The program uses the simple propagator for the Schrodinger equation,
% for a nonlinear time dependent Hamiltonian, of a known functional form.
% For this purpose, it uses a first order approximation, for the propagation to the next time
% step. Hence, It has a limited accuracy.
% The program is intended to a TLS, within the RWA.
% E is a vector of the eigenenergies of the unpurturbed system.
    dt = T/Nt;
    dim = length(chiT);
    chi = zeros(dim, Nt + 1);
    chi(:, Nt + 1) = chiT;
%    field = zeros(1, Nt + 1);
    for ti = (Nt+1):-1:2
        field = 1i/Epenal*(conj(chi(1, ti))*psi(2, ti) - conj(psi(1, ti))*chi(2, ti));
        Ht = [E(1)         -conj(field);
              -field   E(2)            ];
        [P, D] = eig(Ht);
        chit = P\chi(:, ti);
        eigval = diag(D);
        chi(:, ti - 1) = P*(exp(1i*eigval*dt).*chit);
    end
%    field(Nt + 1) = 1i/Epenal*(conj(chi(1, Nt+1))*psi(2, Nt+1) - conj(psi(1, Nt+1))*chi(2, Nt+1));
end