K = input('Enter the number of sinusoidal terms in the line list: ');
dt = input('Enter the time step: ');
N = 2*K;
% an arbitrary initialization of the parameters, d and w:
w = sort(rand(1, K)*2*pi/dt - rand(1, K)*0.5*i);
d = rand(1, K)*10;
% optional - with no damping:
%w = sort(rand(1, K)*3);
c = zeros(N, 1);
for n = 0:(N - 1)
    % The time index is n+1.
    c(n+1) = sum(d.*exp(-i*w*dt*n));
end
% Creating a matrix of the time indices of the signal, c, for U0.
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;
U0 = fft2(c(U0cti));
U1 = fft2(c(U0cti + 1));
[B, Du] = eig(U1, U0);
w
w_result = -log(diag(Du).')/(i*dt);
w_result(real(w_result)<0) = w_result(real(w_result)<0) + 2*pi/dt;
[w_result, orderw] = sort(w_result);
w_result
rwerr = real(w_result - w)./real(w)
imwerr = imag(w_result - w)./imag(w)
werr = sqrt(rwerr.^2 + imwerr.^2)
FTc = fft(c(1:K));
d_result = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
    B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
    d_result(k) = (B(:, k).'*FTc(1:K))^2;
end
d
d_result = d_result(orderw)
derr = real(d_result - d)./d
max_erw = max(abs(rwerr))
max_eimw = max(abs(imwerr))
max_ed = max(abs(derr))
