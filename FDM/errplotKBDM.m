K = input('Enter the number of sinusoidal terms in the line list: ');
N = 2*K;
% an arbitrary initialization of the parameters, d and w:
d = sort(rand(1, K)*10)
w = sort(rand(1, K)*3 - rand(1, K)*3*i)
c = zeros(N, 1);
% Creating a matrix of the time indices of the signal, c, for U0.
v = ones(1, K);
u = 0:K-1;
U0cti = v'*u + u'*v + 1;
allmrwerr = zeros(1, 200);
allmimwerr = zeros(1, 200);
allmderr = zeros(1, 200);
for j = 1:200
    dt = 0.05*j;
    KBDMexi
    allmrwerr(j) = max_erw;
    allmimwerr(j) = max_eimw;
    allmderr(j) = max_ed;
end
plot(0.05:0.05:10, allmrwerr, 'b')
figure
plot(0.05:0.05:10, allmimwerr, 'r')
figure
plot(0.05:0.05:10, allmderr, 'g')
