function U = TDHexpHc(V, ui, xdomain, evdomain, Ncheb, T, Nt)
% The program uses the simple propagator for the Schrodinger equation,
% for a nonlinear time dependent Hamiltonian, of a known functional form.
% For this purpose, it uses a first order approximation, for the propagation to the next time
% step. Hence, It has a limited accuracy.
% V is the potential energy. It's a function handle of the form: @(x, t).
% xdomain = [minx maxx].
% evdomain: The domain of the eigenvalues - [minev maxev].
% Ncheb: The number of Chebychev coefficients in the expansion.
% The program uses the Chebychev expansion for approximation of a function
% of matrix.
    Nx = length(ui);
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
    dt = T/Nt;
    U = zeros(Nx, Nt + 1);
    U(:, 1) = ui;
    leftb = evdomain(1);
    rightb = evdomain(2);
    Ccheb = chebc(@(H) exp(-1i*H*dt), leftb, rightb, Ncheb).';
    for ti = 1:Nt
        U(:, ti + 1) = expHc_1ts(K, V(x, ti*dt), U(:, ti), leftb, rightb, Ncheb, Ccheb);
    end
end