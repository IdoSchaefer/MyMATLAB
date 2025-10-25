% This file assumes the existence of the variables: c (the signal, from t=0), dt (the time step),
% and K (the number of cosine terms in the line list).
% % To find the maximal index needed to the reflection of the signal around 0:
% if (3 - 2*K)<(1 - K)
%     ref_maxi = 2
% Reflection of c around 0:
Nt = length(c);
even_c = [c(Nt:-1:2); c];
% Building the c index matrices:
v = ones(1, K);
u = 0:K-1;
plusi = v.'*u + u.'*v;
minusi = v.'*u - u.'*v;
U0 = even_c(plusi + Nt) + even_c(minusi + Nt);
U1 = 0.5*(even_c(1 + plusi + Nt) + even_c(1 - plusi + Nt) + even_c(1 + minusi + Nt) + even_c(1 - minusi + Nt));
[B, Du] = eig(U1, U0);
w = acos(diag(Du).')/dt;
[temp, w_order] = sort(w);
w = w(w_order);
%w_result = abs(real(expw(1:2:K))) + i*imag(expw(1:2:K))
%w_result = expw(real(expw)>=0);
f = w/(2*pi);
% optional - for grasp results, in cm-1:
wave_number = f*1e15/3e10;
%FTc = fft(c(1:K));
d = zeros(1, K);
for k = 1:K
    % Normalizing B(:, k) with respect to U0:
%    B(:, k) = B(:, k)/sqrt(B(:, k)'*U0*B(:, k));
%    expd(k) = (B(:, k).'*FTc(1:K))^2;
    Bc = (B(:, k).'*c(1:K));
    BU0B = B(:, k)'*U0*B(:, k);
    d(k) = 2*Bc*conj(Bc)*BU0B/(BU0B*conj(BU0B));
end
% wtMat = cos((0:2:(2*K - 1)).'*dt*w);
% d = (wtMat\c(1:2:2*K)).';
d = d(w_order);
% d_result = abs(expd(real(expw)>=0))*2;
% g_result = imag(expw(real(expw)>=0));
tmax = (2*K-1)*dt;
% meand is the mean value of the amplitude, calculated with an integral on
% exp(imag(w)*t) over the time.
% meand = d_result.*1./(g_result*tmax).*(exp(g_result*tmax) - 1);
%d_result = abs(expd(1:2:K))*2
%[maxd, nmaxd] = max(abs(expd)*2);
[maxd, nmaxd] = max(d);
%wn_maxd = abs(real(w_result(nmaxd))*1e15/(3e10*2*pi))
%wreal = real(w_result);
%wn_maxd = real(w_result(nmaxd))*1e15/(3e10*2*pi)
w_maxd = w(nmaxd)
maxd
figure
%plot(wave_number, d_result)
plot(w, d)
hold on
%plot(wave_number, meand, '--r')