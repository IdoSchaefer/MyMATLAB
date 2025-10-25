L = 2*pi;
max_level = 4;
N = 2^max_level;
h = L/N;
omega = 0.5;
x = (h:h:(L - h)).';
n = 5;
fx = sin(n*x/2);
f = n^2/4*kron(fx, fx);
Lmat = Laplacian_matrix(N, h);
Phi_analytic = -2/n^2*f;
Phi_matlab = Lmat\f;
Phi_guess = zeros((N - 1)^2, 1);
Phi_guess(ceil(end/2)) = 1;
Niter = 10;
Ncycles = 20;
normres = zeros(Ncycles + 1, 1);
error = zeros(Ncycles + 1, 1);
error_num = zeros(Ncycles + 1, 1);
normres(1) = norm(f - Lmat*Phi_guess)/(N - 1);
error(1) = norm(Phi_guess - Phi_analytic)/norm(Phi_analytic);
error_num(1) = norm(Phi_guess - Phi_matlab)/norm(Phi_matlab);
Phi_out = Phi_guess;
for cyclei = 1:Ncycles
    Phi_out = multigrid(omega, h, Phi_out, f, Niter, 1);
    normres(cyclei + 1) = norm(f - Lmat*Phi_out)/(N - 1);
    error(cyclei + 1) = norm(Phi_out - Phi_analytic)/norm(Phi_analytic);
    error_num(cyclei + 1) = norm(Phi_out - Phi_matlab)/norm(Phi_matlab);
end
figure
plot(0:Ncycles, log10(normres))
figure
plot(0:Ncycles, log10(error))
figure
plot(0:Ncycles, log10(error_num))