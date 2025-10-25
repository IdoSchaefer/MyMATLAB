for n = 0:(N - 1)
    % The time index is n+1.
    c(n+1) = sum(d.*exp(-i*w*dt*n));
end
U0 = c(U0cti);
U1 = c(U0cti + 1);
[B, Du] = eig(U1, U0);
w_result = (-log(diag(Du).')/(i*dt));
w_result(real(w_result)<0) = w_result(real(w_result)<0) + 2*pi/dt;
[w_result, orderw] = sort(w_result);
rwerr = real(w_result - w)./real(w);
imwerr = imag(w_result - w)./imag(w);
%werr = sqrt(rwerr.^2 + imwerr.^2);
d_result = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
    B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
    d_result(k) = (B(:, k).'*c(1:K))^2;
end
d_result = d_result(orderw);
derr = real(d_result - d)./d;
max_erw = max(abs(rwerr))
max_eimw = max(abs(imwerr))
max_ed = max(abs(derr))
