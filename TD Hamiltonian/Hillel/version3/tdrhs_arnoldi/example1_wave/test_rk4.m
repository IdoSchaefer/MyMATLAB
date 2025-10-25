function test_rk4
global L
global x
global T
global w0
global dx
n=64;
dx=1/n;
for j=1:n
    x(j,1)=1-(j-1)*dx;
end
u=zeros(n,1);
T=1.5;
steps=1000;
dt=T/steps
t=0;
maxit=100;
w0=1;
u0=u;
tic
for k=1:steps
    u=rk4(u,'matvecrk4',t,dt);
    t=t+dt;
    u(1)=sin(t/1);
    u0=u;
end
toc
plot(x,u)
u1=u0;
save fu1 u1
load fu1 u1;
norm(u1-u0)


