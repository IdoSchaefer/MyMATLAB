function [psi, field] = TLSKpsi(E, psi0, chi, Epenal, T, Nt)
% The program uses the simple propagator for the Schrodinger equation,
% for a nonlinear time dependent Hamiltonian, of a known functional form.
% For this purpose, it uses a first order approximation, for the propagation to the next time
% step. Hence, It has a limited accuracy.
% The program is intended to a TLS, within the RWA.
% E is a vector of the eigenenergies of the unpurturbed system.
    dt = T/Nt;
    dim = length(psi0);
    psi = zeros(dim, Nt + 1);
    psi(:, 1) = psi0;
    field = zeros(1, Nt + 1);
    for ti = 1:Nt
        field(ti) = 1i/Epenal*(conj(chi(1, ti))*psi(2, ti) - conj(psi(1, ti))*chi(2, ti));
        Ht = [E(1)         -conj(field(ti));
              -field(ti)   E(2)            ];
        [P, D] = eig(Ht);
        psit = P\psi(:, ti);
        eigval = diag(D);
        psi(:, ti + 1) = P*(exp(-1i*eigval*dt).*psit);
    end
    field(Nt + 1) = 1i/Epenal*(conj(chi(1, Nt+1))*psi(2, Nt+1) - conj(psi(1, Nt+1))*chi(2, Nt+1));
end