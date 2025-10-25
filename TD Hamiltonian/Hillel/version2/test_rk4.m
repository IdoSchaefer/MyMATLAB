function test_rk4
global L
global x
global T
global w0
n=128;
[x,u,L]=initial(n);
T=15;
steps=3000;
dt=T/steps
t=0;
maxit=100;
w0=1;
u0=u;
tic
for k=1:steps
    u=rk4(u,t,dt);
    t=t+dt;
    hold on
    nu=norm(u-u0)
    u0=u;
end
toc
plot(abs(u))
% u1=u0;
% save fu1 u1
load fu1 u1;
norm(u1-u0)


