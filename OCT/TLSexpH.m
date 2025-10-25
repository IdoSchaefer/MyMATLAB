function [U, field] = TLSexpH(E, fieldfun, ui, T, Nt)
% The program uses the simple propagator for the Schrodinger equation,
% for a nonlinear time dependent Hamiltonian, of a known functional form.
% For this purpose, it uses a first order approximation, for the propagation to the next time
% step. Hence, It has a limited accuracy.
% The program is intended to a TLS, within the RWA.
% E is a vector of the eigenenergies of the unpurturbed system.
% fieldfun represents the functional form of the field. 
% It's a function handle of the form: @(u, t).
    dt = T/Nt;
    dim = length(ui);
    U = zeros(dim, Nt + 1);
    U(:, 1) = ui;
    field = zeros(1, Nt + 1);
    for ti = 1:Nt
        field(ti) = fieldfun(U(:, ti), (ti-1)*dt);
        Ht = [E(1)         -conj(field(ti));
              -field(ti)   E(2)            ];
        [P, D] = eig(Ht);
        ut = P\U(:, ti);
        eigval = diag(D);
        U(:, ti + 1) = P*(exp(-1i*eigval*dt).*ut);
    end
    field(Nt + 1) = fieldfun(U(:, Nt + 1), T);
end