K = input('Enter the number of sinusoidal terms in the line list: ');
dt = input('Enter the time step: ');
N = 2*K;
M = K;
% an arbitrary initialization of the parameters, d and w:
d = sort(rand(1, K)*10);
w = sort(rand(1, K)*2*pi/dt - rand(1, K)*0.5*i);
% optional - with no damping:
%w = sort(rand(1, K)*pi/dt);
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
Kb = 10;
stepb = 6;
nuse = ceil(Kb/2) - ceil(stepb/2) + 1;
w_result = zeros(1, M);
FTc = fft(c(1:K));
d_result = zeros(1, M);
jj = 1;
getwb
% sorting the data to get the lowest frequencies in the beginning of the block:
wsort
w_result(1:(nuse + stepb - 1)) = wb(1:(nuse + stepb - 1));
for k = 1:nuse + stepb - 1
    % Normalizing Bb(:, k) with respect to U0:
    Bb(:, k) = Bb(:, k)/sqrt(Bb(:, k).'*U0b*Bb(:, k));
    d_result(k) = (Bb(:, k).'*FTc(1:Kb))^2;
end
for jj = 1 + stepb:stepb:M - Kb + 1
    getwb
% sorting the data to get the middle frequencies in the middle of the block:
    wsort
    w_result((jj + nuse - 1):(jj + nuse + stepb - 2)) = wb(nuse:nuse + stepb - 1);
    for k = nuse:nuse + stepb - 1
        % Normalizing Bb(:, k) with respect to U0b:
        Bb(:, k) = Bb(:, k)/sqrt(Bb(:, k).'*U0b*Bb(:, k));
        d_result(jj + k - 1) = (Bb(:, k).'*FTc(jj:Kb + jj - 1))^2;
    end
end
leftw = rem(M - Kb, stepb); 
if leftw > 0
    jj = M - Kb + 1;
    getwb
    wsort
%    w_result(M - Kb + nuse:M) = wb(nuse:Kb)
end
    w_result(M - Kb + nuse + stepb - leftw:M) = wb(nuse + stepb - leftw:Kb)

for k = nuse + stepb - leftw:Kb
    % Normalizing Bb(:, k) with respect to U0b:
    Bb(:, k) = Bb(:, k)/sqrt(Bb(:, k).'*U0b*Bb(:, k));
    d_result(M - Kb + k) = (Bb(:, k).'*FTc(M - Kb + 1:M))^2;
end
w
%w_result(real(w_result)<0) = w_result(real(w_result)<0) + 2*pi/dt;
w_result = sort(w_result)
rwerr = real(w_result - w)./real(w)
imwerr = imag(w_result - w)./imag(w)
werr = sqrt(rwerr.^2 + imwerr.^2)
d
d_result = sort(d_result)
derr = real(d_result - d)./d
max_erw = max(abs(rwerr))
max_eimw = max(abs(imwerr))
max_ed = max(abs(derr))
