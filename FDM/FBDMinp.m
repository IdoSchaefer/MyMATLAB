% This file assumes the existance of the variables: c (the signal, from t=0), dt (the time step),
% and K (the number of complex exponential terms in the line list).
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;
U0 = fft2(c(U0cti));
U1 = fft2(c(U0cti + 1));
[B, Du] = eig(U1, U0);
w_result = -log(diag(Du).')/(i*dt);
w_result(real(w_result)<0) = w_result(real(w_result)<0) + 2*pi/dt;
[temp, rw_order] = sort(real(w_result));
w_result = w_result(rw_order)
%w_result = abs(real(expw(1:2:K))) + i*imag(expw(1:2:K))
f = real(w_result)/(2*pi);
% optional - for grasp results, in cm-1:
wave_number = f*1e15/3e10
FTc = fft(c(1:K));
expd = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
    B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
    expd(k) = (B(:, k).'*FTc(1:K))^2;
end
expd = expd(rw_order)
d_result = abs(expd)
[maxd, nmaxd] = max(d_result);
%wn_maxd = abs(real(w_result(nmaxd))*1e15/(3e10*2*pi))
wn_maxd = real(w_result(nmaxd))*1e15/(3e10*2*pi)
maxd