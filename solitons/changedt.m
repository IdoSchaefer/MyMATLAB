dt = input('Enter a new dt: ');
t = 0:dt:tfinal;
x = deval(sol, t);
toptimei = fix(tfinal/dt + 1);
timei = 1:toptimei;
allrvalues = x(2:N, timei) - x(1:(N-1), timei);
maxr = max(max(allrvalues));
minr = min(min(allrvalues));