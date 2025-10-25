% Constructing the modified Chebyshev grid:
L = 16;
N = 256;
theta = (0:N).'*pi/N;
xcheb =  cos(theta);
%p = pi/2*0.955;
p = 2*atan((1e-3)^(1/(N + 1)));
x = 1/p*asin(xcheb*sin(p))*0.5*L;
% Weight function for integration:
w = sin(p)/p*L/2*sin(theta)./sqrt(1 - xcheb.^2*sin(p)^2)*pi/N;
% x = flipud(x);
T = 10;
Nt = 100;
dt = T/Nt;
tpoints = 0:dt:T;
% Dxinv is 1/g'(y; p):
% Dxinv = p*sqrt(1-(sin(p)*xcheb).^2)/sin(p);
% Or, alternatively, dy/dx*L/2:
% Dxinv = (p/sin(p))*cos(2*p*x/L);
% The initial state (unnormalized):
psi0 = exp(1i*x).*pi^(-1/4).*exp(-x.^2/2);
psi0 = psi0/sqrt(psi0'*(w.*psi0));
% The potential:
V = x.^2/2;
% The coefficients for the 2'nd derivative:
[D2coef, Dcoef] = D2modchebbcoef(N, L, p);
% The Hamiltonian operation:
Hop = @(psi) Hpsimodchebp(V, psi, L, D2coef, Dcoef, 1);
% Finding the eigenvalue domain:
H = Hmodchebb(@(x) x.^2/2, L, N, p, 1);
[leftb, rightb, ev] = Re_evdomain(H);
% Propagating for a short time interval - 0.1, using the Chebyshev algorithm:
%Psi = fMchebop(Hop, leftb, rightb, @(H) exp(-1i*H), psi0, T, Nt, 1e4);
Psi = fMchebop_tol(Hop, leftb, rightb, @(H) exp(-1i*H), psi0, T, Nt, 1e4, 1e-10);
figure
plot(x, conj(Psi(:, end)).*Psi(:, end))
mx = evmiugen(Psi, x, 1, w);
er = mx - sin(tpoints);
figure
plot(tpoints, er)
%er = sin(10) - (Psi(:, end)'*(x.*w.*Psi(:, end))) %/(Psi(:, end)'*(w.*Psi(:, end)))
% Propagating for a short time interval - 0.1, using the simple Arnoldi 
% algorithm:
% PsiKr = fMkropsimp(Hop, @(H) exp(-1i*H), psi0, 0.1, 100, 100);
% figure
% % The result diverges again, but at a slower rate:
% plot(x, conj(PsiKr(:, end)).*PsiKr(:, end))
% erKr = sin(0.1) - (PsiKr(:, end)'*(x.*w.*PsiKr(:, end)))/(PsiKr(:, end)'*(w.*PsiKr(:, end)))
