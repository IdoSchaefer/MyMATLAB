if boundary == 'f'
    sol = ode45(@fixedchain, [0 tfinal], [initx initv], [], F, r, m, N);
elseif boundary == 'c'
    sol = ode45(@cyclicchain, [0 tfinal], [initx initv], [], F, r, r0, m, N);
else
    sol = ode45(@openchain, [0 tfinal], [initx initv], [], F, r, m, N);
end
dt = 0.01;
t = 0:dt:tfinal;
x = deval(sol, t);
toptimei = fix(tfinal/dt + 1);
timei = 1:toptimei;
allrvalues = x(2:N, timei) - x(1:(N-1), timei);
maxr = max(max(allrvalues));
minr = min(min(allrvalues));
tstep = 0.05;
tistep = 5;
clf
hold on
i = 1:N;
plot(t, x(i,:))