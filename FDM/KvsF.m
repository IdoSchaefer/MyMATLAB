K = input('Enter the number of sinusoidal terms in the line list: ');
dt = input('Enter the time step: ');
N = 2*K;
Ncheck = 50;
errdiff = zeros(3, Ncheck);
for checki = 1:Ncheck
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
    U0 = c(U0cti);
    U1 = c(U0cti + 1);
    U2w
    U2d
    max_erwK = max(abs(rwerr));
    max_eimwK = max(abs(imwerr));
    max_edK = max(abs(derr));
    U0 = fft2(c(U0cti));
    U1 = fft2(c(U0cti + 1));
    U2w
    c = fft(c(1:K));
    U2d    
    max_erwF = max(abs(rwerr));
    max_eimwF = max(abs(imwerr));
    max_edF = max(abs(derr));
    errdiff(1, checki) = max_erwK - max_erwF;
    errdiff(2, checki) = max_eimwK - max_eimwF;
    errdiff(3, checki) = max_edK - max_edF;
end
plot(1:Ncheck, errdiff(1, :))
figure
plot(1:Ncheck, errdiff(2, :))
figure
plot(1:Ncheck, errdiff(3, :))