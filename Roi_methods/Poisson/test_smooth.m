L = 2*pi;
%N = 8;
N = 16;
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
Niter = 100;
allnormr = zeros(5, Niter + 1);
error = zeros(5, 1);
error_num = zeros(5, 1);
for omegai = 1:5
    omega = (omegai - 1)*0.5;
    [Phi_out, allnormr(omegai, :)] = smooth(update_matrix(N, omega), Lmat, omega, h, Phi_guess, f, Niter);
    error(omegai) = norm(Phi_out - Phi_analytic)/norm(Phi_analytic);
    error_num(omegai) = norm(Phi_out - Phi_matlab)/norm(Phi_matlab);
end
figure
plot(0:Niter, allnormr)
figure
plot((0:4)*0.5, error)
figure
plot((0:4)*0.5, error_num)