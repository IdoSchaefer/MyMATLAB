K = input('Enter the number of sinusoidal terms in the line list (K): ');
sdt = input('Enter the time step for N=2*K: ');
% an arbitrary initialization of the parameters, d and w:
d = sort(rand(1, K)*10);
w = sort(rand(1, K)*3 - rand(1, K)*3*i);
maxN = 20*K;
allmed = zeros(1, (maxN - (2*K))/2 + 1);
for j = 1:(maxN-(2*K))/2 + 1
    N = 2*K + 2*j - 2;
    M = N/2;
    dt = 2*K*sdt/N;
    c = zeros(N, 1);
    for n = 0:(N - 1)
        % The time index is n+1.
        c(n+1) = sum(d.*exp(-i*w*dt*n));
    end
% Creating a matrix of the time indices of the signal, c, for U0.
    v = ones(1, M);
    u = 0:M-1;
    U0cti = v'*u + u'*v + 1;
    U0 = c(U0cti);
    U1 = c(U0cti + 1);
    [B, Du] = eig(U1, U0);
%w(M-K+1:M) = w;
%w(1:M-K) = 0
    %w_result = sort(-log(diag(Du).')/(i*dt))
%rwerr = real(w_result - w)./real(w)
%imwerr = imag(w_result - w)./imag(w)
%werr = sqrt(rwerr.^2 + imwerr.^2)
    d_result = zeros(1, M);
    for k = 1:M
    % Normalizing B(:, k) with respect to U0:
        B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
        d_result(k) = (B(:, k).'*c(1:M))^2;
    end
%d(M-K+1:M) = d;
%d(1:M-K) = 0
%d
    d_result = sort(d_result);
%d_result(M-K+1:M)
    derr = real(d_result(M-K+1:M) - d)./d;
%max_erw = max(abs(rwerr))
%max_eimw = max(abs(imwerr))
    max_ed = max(abs(derr));
    allmed(j) = max_ed;
end
plot((1:(maxN-(2*K))/2 + 1)*2 + 2*K - 2, allmed)