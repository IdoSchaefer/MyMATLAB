% Constructing the Chebyshev grid:
L = 16;
N = 128;
theta = ((1:N).'*2 - 1)*pi/(2*N);
xcheb = cos(theta);
x = 0.5*xcheb*L;
% The initial state (unnormalized):
psi0 = exp(1i*x).*pi^(-1/4).*exp(-x.^2/2);
% The potential:
V = x.^2/2;
% The Hamiltonian operation:
Hop = @(psi) Hpsichebp(V, psi, L, 1);
% Propagating for a short time interval - 1, using the Chebyshev algorithm:
Psi = fMchebop(Hop, -340, 100, @(H) exp(-1i*H), psi0, 1, 100, 1000);
figure
% The result diverges:
plot(x, Psi(:, end).*conj(Psi(:, end))/(Psi(:, end)'*(Psi(:, end).*sin(theta))));
% Propagating for a short time interval - 0.1, using the simple Arnoldi 
% algorithm:
PsiKr = fMkropsimp(Hop, @(H) exp(-1i*H), psi0, 0.1, 100, 100);
figure
% The result diverges again, but at a faster rate:
plot(x, PsiKr(:, end))

