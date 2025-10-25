function U = TDHexpH(E, Vt, ui, T, Nt)
% The program uses the simple propagator for the Schrodinger equation,
% for a nonlinear time dependent Hamiltonian, of a known functional form.
% For this purpose, it uses a first order approximation, for the propagation to the next time
% step. Hence, It has a limited accuracy.
% The program is intended to a TLS, within the RWA.
% E is a vector of the eigenenergies of the unpurturbed system.
% Vt represent the purturbation. It is a 3d tensor, when the 3'rd index
% represent the time index, and the first 2 indices are the matrix. 
    dt = T/Nt;
    dim = length(ui);
    U = zeros(dim, Nt + 1);
    U(:, 1) = ui;
    for ti = 1:Nt
        Ht = diag(E) + Vt(:, :, ti);
        [P, D] = eig(Ht);
        ut = P\U(:, ti);
        eigval = diag(D);
        U(:, ti + 1) = P*(exp(-1i*eigval*dt).*ut);
    end
end