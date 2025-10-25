Nt = 100;
K = Nt/2;
dt = 0.5;
c = zeros(Nt, 1);
%dtest = (rand(1, Kre) - 0.5)*5;
dtest = (rand(1, K) - 0.5)*5;
wtest = sort(rand(1, K)*pi/dt);
for n = 0:(Nt - 1)
    % The time index is n+1.
    c(n+1) = sum(dtest.*cos(wtest*dt*n));
end
KBDMr
%der = d(abs(d)>min(abs(dtest))) - dtest(dtest>min(abs(dtest)));
der = d - dtest;
wer = w - wtest;
max_ed = max(abs(der))
max_ew = max(abs(wer))
