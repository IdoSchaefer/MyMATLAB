function mE = evH(psi, Hpsif)
% The function computes the time independent Hamiltonian expectation value in each time step
% of psi.
% Hpsif: a function handle; returns the operation of H on psi.
    Nt = size(psi, 2);
    mE = zeros(1, Nt);
    for ti = 1:Nt
        mE(ti) = psi(:, ti)'*Hpsif(psi(:, ti));
    end
end