d_result = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
    B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
    d_result(k) = (B(:, k).'*c(1:K))^2;
end
d_result = d_result(orderw);
derr = real(d_result - d)./d;
