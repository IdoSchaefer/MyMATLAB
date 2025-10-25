function test
global L
global x
global T
global w0
global matvecs
matvecs=0;
n=128;
[x,u,L]=initial(n);
m1=9;
m2=9;
steps=300;
T=15;
dt=T/steps;
t=0;
w0=1;
z1=0;z2=-sqrt(-1)*250;
tic
for k=1:steps
    u=tdrhs_chb(u,'matvec',t,dt,m1,m2,z1,z2);
    t=t+dt;
end
toc
plot(abs(u(:,1)))
load fu1 u1
norm(u1-u(:,1))
matvecs

