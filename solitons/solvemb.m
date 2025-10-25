clf
syms F V r;
funchoice = input('Choose: a potential/force function(V/F):', 's');
if funchoice == 'V'
        getV = input('Enter a potential function in the terms of r: ', 's');
        V = eval(getV);
        F = -diff(V);
    else
        getF = input('Enter a force function in the terms of r: ', 's');
        F = eval(getF);
end
N = input('Enter the number of the bodies: ');
m = input('Enter the mass:');
tfinal = input('Enter the time interval:');
solveF = subs(solve(F));
r0 = solveF(1);
initx = 0:r0:(N-1)*r0;
initv = zeros(1, N);
boundary = input('Choose boundary conditions: open, fixed, or cyclic (o, f, c):', 's');
initchoice = input('Choose: initialize values as your choice, as a sinus or as a gaussian(c, s, g)?', 's');
    if initchoice ~= 's' && initchoice ~= 'g'
        num = input('Enter the number of the body:');
        while num>=1 && num<=N
            initx(num) = input('Enter the initial x:');
            initv(num) = input('Enter the initial v:');
            num = input('Enter the number of the body:');
        end
    elseif initchoice == 'g'
        A = input('Enter the gaussians height: ');
        sigma = input('Enter the standard deviation, with the body number as the variable:');
        n0 = input('Enter the center of the gaussian, with the body number as the variable:');
        ri = 1:(N-1);
        rvalue = r0 - A*exp(-(ri - n0).^2/(2*sigma^2));
        rfit = ((N-1)*r0 - sum(rvalue))/(N-1);
        rvalue = rvalue + rfit;
        for xi = 2:N
            initx(xi) = initx(xi-1) + rvalue(xi - 1);
        end
    end
if boundary == 'f'
    if initchoice == 's'   
        mode = input('Which mode?');
        A = input('Enter the amplitude:');
        q = A*sin(pi*mode/((N-1)*r0)*(0:r0:(N-1)*r0));
        initx = initx + q;
    end
    sol = ode45(@manybodies1, [0 tfinal], [initx initv], [], F, r, m, N);
elseif boundary == 'c'
    if initchoice == 's'
        Ncycles = input('How many cycles?');
        A = input('Enter the amplitude:');
        q = A*sin(2*pi*Ncycles/(N*r0)*(0:r0:(N-1)*r0));
        initx = initx + q;
    end
    sol = ode45(@cyclicmanybodies, [0 tfinal], [initx initv], [], F, r, r0, m, N);
else
    if initchoice == 's'
        mode = input('Which mode?');
        A = input('Enter the amplitude:');
        q = A*sin(pi*mode/((N-1)*r0)*(0:r0:(N-1)*r0));
        initx = initx + q;
    end
    sol = ode45(@manybodies, [0 tfinal], [initx initv], [], F, r, m, N);
end
t = 0:0.1:tfinal;
x = deval(sol, t);
hold on
i = 1:N;
plot(t, x(i,:))
figure
axtop = max(x(N,:));
axbot = min(x(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2]);
if N<=32
    ballsize = 24;
else
    ballsize = (32/N)*24;
end
for j = 1:N
    ball(j) = line(0, 0, 'marker', '.', 'markersize', ballsize);
end
toptimei = fix(tfinal/0.1 + 1);
for timei = 1:toptimei
    for j=1:N
        set(ball(j), 'xdata', x(j, timei))
    end
    pause(0.1)
end