clf
syms F V r;
getchoice = input('Choose: a potential/force function(V/F):', 's');
if getchoice == 'V'
        getV = input('Enter a potential function in the therms of r: ', 's');
        V = eval(getV);
        F = -diff(V);
    else
        getF = input('Enter a force function in the therms of r: ', 's');
        F = eval(getF);
end
N = input('Enter the number of the bodies: ');
m = input('Enter the mass:');
tfinal = input('Enter the time interval:');
r0 = subs(solve(F));
initx = 0:r0:(N-1)*r0;
initv = zeros(1, N);
getchoice = input('Choose: initialize values as your choice, or as a sinus (c,s)?', 's');
if getchoice ~= 's'
        num = input('Enter the number of the body (the outer bodies are fixed):');
        while num>=2 && num<=(N-1)
            initx(num) = input('Enter the initial x:');
            initv(num) = input('Enter the initial v:');
            num = input('Enter the number of the body (the outer bodies are fixed):');
        end
    else
        mode = input('Which mode?');
        A = input('Enter the amplitude:');
        q = A*sin(pi*mode/((N-1)*r0)*(0:r0:(N-1)*r0));
        initx = initx + q;
end
xv = frog(@manybodiesf, [initx initv], [0 tfinal], 1e-3);
hold on
i = 1:N;
plot(t, xv(i,:))
figure
axtop = max(xv(N,:));
axbot = min(xv(1,:));
axis([axbot - 1 , axtop + 1, -(axtop - axbot)/2, (axtop - axbot)/2]);
for j = 1:N
    ball(j) = line(0, 0, 'marker', '.', 'markersize', 24);
end
    toptimei = fix(tfinal/0.1 + 1);
for timei = 1:toptimei
    for j=1:N
        set(ball(j), 'xdata', x(j, timei))
    end
    pause(0.1)
end

%Nested function.
function DxDv = manybodiesf(getxv)
    Dx = zeros(N, 1);
    Dv = zeros(N, 1);
    ibody = 1:N;
    iforce = 2:N;
    x = getxv(ibody);
    v = getxv(N + ibody);
    Dx(ibody) = v(ibody);
    rvalue(iforce) = x(iforce) - x(iforce-1);
    Fvalue = zeros(1, N+1);
    Fvalue(iforce) = subs(F, r, rvalue(iforce));
    Fvalue(N+1) = Fvalue(2);
    Dv(1) = 0;
    ibody = 2:(N-1);
    Dv(ibody) = (-Fvalue(ibody+1)+Fvalue(ibody))/m;
    Dv(N) = 0;
    DxDv = [Dx; Dv];
end

end