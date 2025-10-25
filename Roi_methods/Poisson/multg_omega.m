L = 2*pi;
max_level = 4;
N = 2^max_level;
h = L/N;
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
%allnormr = zeros(5, Niter + 1);
normres = zeros(5, 1);
error = zeros(5, 1);
error_num = zeros(5, 1);
for omegai = 1:5
    omega = (omegai - 1)*0.5;
    Phi_out = Phi_guess;
    for cyclei = 1:Ncycles
        Phi_out = multigrid(omega, h, Phi_out, f, Niter, 1);
        %Phi_out = smooth(update_matrix(N, omega), Lmat, omega, h, Phi_out, f, Niter);
    end
    normres(omegai) = norm(f - Lmat*Phi_out)/(N - 1);
    error(omegai) = norm(Phi_out - Phi_analytic)/norm(Phi_analytic);
    error_num(omegai) = norm(Phi_out - Phi_matlab)/norm(Phi_matlab);
end
% figure
% plot(0:Niter, allnormr)
figure
plot((0:4)*0.5, log10(normres))
figure
plot((0:4)*0.5, log10(error))
figure
plot((0:4)*0.5, log10(error_num))