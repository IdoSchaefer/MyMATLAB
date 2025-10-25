K = input('Enter the number of sinusoidal terms in the line list: ');
dt = input('Enter the time step: ');
N = 2*K;
% an arbitrary initialization of the parameters, d and w:
d = sort(rand(1, K)*10);
w = sort(rand(1, K)*2*pi/dt - rand(1, K)*0.5*i);
c = zeros(N, 1);
% Creating a matrix of the time indices of the signal, c, for U0.
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;