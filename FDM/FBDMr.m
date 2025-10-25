Kre = input('Enter the number of sinusoidal terms in the line list: ');
% Each (real) cos function is spanned by 2 complex sinusoidal terms:
K = 2*Kre;
dt = input('Enter the time step: ');
N = 2*K;
% an arbitrary initialization of the parameters, d and w:
d = rand(1, Kre)*10;
% Creating a real signal.
rw = sort(rand(1, Kre)*pi/dt);
g = rand(1, Kre)*0.5;
c = zeros(N, 1);
for n = 0:(N - 1)
    % The time index is n+1.
    c(n+1) = sum(d.*cos(rw*dt*n).*exp(-g*dt*n));
end
% Creating a matrix of the time indices of the signal, c, for U0.
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;
U0 = fft2(c(U0cti));
U1 = fft2(c(U0cti + 1));
[B, Du] = eig(U1, U0);
w = rw - i*g
expw = -log(diag(Du).')/(i*dt);
[temp, rw_order] = sort(abs(real(expw)));
expw = expw(rw_order)
%w_result = abs(real(expw(1:2:K))) + i*imag(expw(1:2:K))
w_result = expw(real(expw)>=0)
rwerr = real(w_result - rw)./rw
imwerr = imag(w_result - w)./imag(w)
werr = sqrt(rwerr.^2 + imwerr.^2)
FTc = fft(c(1:K));
expd = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
    B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
    expd(k) = (B(:, k).'*FTc(1:K))^2;
end
d
expd = expd(rw_order)
d_result = abs(real((expd(1:2:K)))*2
derr = real(d_result - d)./d
max_erw = max(abs(rwerr))
max_eimw = max(abs(imwerr))
max_ed = max(abs(derr))
