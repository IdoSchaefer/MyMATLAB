clf
syms F V r
getchoice = input('Choose: a potential/force function(V/F):', 's');
if getchoice == 'V'
        getV = input('Enter a potential function in the therms of r: ', 's');
        V = eval(getV);
        F = -diff(V);
    else
        getF = input('Enter a force function in the therms of r: ', 's');
        F = eval(getF);
end
m = input('Enter the mass:');
tfinal = input('Enter the time interval:');
r0 = subs(solve(F));
initx = [0 r0];
initv = zeros(1, 2);
num = input('Enter the number of the body:');
while num>=1 && num<=2
    initx(num) = input('Enter the initial x:');
    initv(num) = input('Enter the initial v:');
    num = input('Enter the number of the body:');
end
sol = ode45(@twobodies, [0 tfinal], [initx initv], [], F, r, m);
t = 0:0.1:tfinal;
x = deval(sol, t);
plot(t, x(1,:), 'b')
hold on
plot(t, x(2,:), 'k')
figure
axtop = max(x(2,:));
axbot = min(x(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2]);
ball(1) = line(0, 0, 'marker', '.', 'markersize', 24);
ball(2) = line(0, 0, 'marker', '.', 'markersize', 24);
topi = fix(tfinal/0.1 + 1);
for i=1:topi
    set(ball(1), 'xdata', x(1, i))
    set(ball(2), 'xdata', x(2, i))
    pause(0.1)
end