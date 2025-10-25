Nt = 100;
K = Nt/2;
dt = 1;
c = zeros(Nt, 1);
%dtest = (rand(1, Kre) - 0.5)*5;
dtest = exp(-1i*rand(1, K)*2*pi);
wtest = sort((rand(1, K) - 0.5)*2*pi/dt);
for n = -(K-1):K
    % The time index is n+1.
    c(n+K) = sum(dtest.*cos(wtest*dt*n));
end
KBDMrw
%der = d(abs(d)>min(abs(dtest))) - dtest(dtest>min(abs(dtest)));
der = d - dtest;
wer = w - wtest;
max_ed = max(abs(der))
max_ew = max(abs(wer))