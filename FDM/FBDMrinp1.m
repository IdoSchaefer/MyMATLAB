% This file assumes the existance of the variables: c (the signal, from t=0), dt (the time step),
% and Kre (the number of real sinusoidal terms in the line list).
K = 2*Kre;
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;
allc = c;
maxci = length(c) - 2*K + 1;
alld = zeros(1, maxci);
allwn = zeros(1, maxci);
for ci = 1:maxci
    c = allc(ci:(ci + 2*K - 1));
    U0 = fft2(c(U0cti));
    U1 = fft2(c(U0cti + 1));
    [B, Du] = eig(U1, U0);
    expw = -log(diag(Du).')/(i*dt);
%[temp, rw_order] = sort(abs(real(expw)));
%expw = expw(rw_order)
%w_result = abs(real(expw(1:2:K))) + i*imag(expw(1:2:K))
%w_result = expw(real(expw)>=0)
%f = real(w_result)/(2*pi);
% optional - for grasp results, in cm-1:
%wave_number = f*1e15/3e10
    FTc = fft(c(1:K));
    expd = zeros(1, K);
    for k = 1:K
        % Normalizing B(:, k) with respect to U0:
        B(:, k) = B(:, k)/sqrt(B(:, k).'*U0*B(:, k));
        expd(k) = (B(:, k).'*FTc(1:K))^2;
    end
%expd = expd(rw_order)
%d_result = abs(expd(1:2:K))*2
    [maxd, nmaxd] = max(abs(expd)*2);
    wn_maxd = abs(real(expw(nmaxd))*1e15/(3e10*2*pi));
    alld(ci) = maxd;
    allwn(ci) = wn_maxd;
end
meanwn = mean(allwn)
sigma = std(allwn)
maxwn = max(allwn)
minwn = min(allwn)
% plot(1:maxci, allwn)
% optional, instead of last row:
plot((0:maxci-1)*dt, allwn)