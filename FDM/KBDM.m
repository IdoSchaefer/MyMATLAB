K = input('Enter the number of sinusoidal terms in the line list: ');
N = 2*K;
% an arbitrary initialization of the parameters, d and w:
d = 1:K;
w = (1:K)*(1 - i);
c = zeros(N, 1);
for n = 0:(N - 1)
    % The time index is n+1.
    % tau = 1
    c(n+1) = sum(d.*exp(-i*w*n));
end
% Creating a matrix of the time indices of the signal, c, for U0.
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;
U0 = c(U0cti);
U1 = c(U0cti + 1);
[B, Du] = eig(U1, U0);
w_result = -log(diag(Du).')/i
d_result = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
    B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
    d_result(k) = (B(:, k).'*c(1:K))^2;
end
d_result
