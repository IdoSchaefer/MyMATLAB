%%% Solves newton's second law for an x dependent force function.
figure
F = input('Enter a force function in the therms of x and v:', 's');
m = input('Enter mass:');
tfinal = input('Enter the time range:');
initx = input('Enter the initial x position:');
initv = input('Enter the initial velocity:');
sol = ode45(@newton, [0 tfinal], [initx initv], [], F, m);
t = 0:0.1:tfinal;
x = deval(sol, t);
plot(t, x(1,:))
figure
axtop = max(x(1,:));
axbot = min(x(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2])
ball=line(0, 0, 'marker', '.', 'markersize', 24);
topi = fix(tfinal/0.1 + 1);
for i=1:topi
    set(ball, 'xdata', x(1, i))
    pause(0.1)
end