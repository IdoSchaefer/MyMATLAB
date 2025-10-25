function U = TDHexpH1(V, ui, xdomain, T, Nt)
%%% The program isn't effective - for small matrices of many level systems,
%%% use TDHexpH. For large matrices, use TDHexpHc.
% The program uses the simple propagator for the Schrodinger equation,
% for a nonlinear time dependent Hamiltonian, of a known functional form.
% For this purpose, it uses a first order approximation, for the propagation to the next time
% step. Hence, It has a limited accuracy.
% V is the potential energy. It's a function handle of the form: @(x, t).
    Nx = length(ui);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    dp = 2*pi/xdlength;
    p = (-pi/dx):dp:(pi/dx - dp);
% The kinetic energy matrix in the p domain:
    Kp = diag(p.^2/2);
% The kinetic energy matrix in the x domain:
    K = 128*ifft(ifft(ifftshift(Kp))')';
    dt = T/Nt;
    U = zeros(Nx, Nt + 1);
    U(:, 1) = ui;
    for ti = 1:Nt
% The Hamiltonian:
        Ht = diag(V(x, ti*dt)) + K;
        [P, D] = eig(Ht);
        ut = P\U(:, ti);
        eigval = diag(D);
        U(:, ti + 1) = P*(exp(-1i*eigval*dt).*ut);
    end
end