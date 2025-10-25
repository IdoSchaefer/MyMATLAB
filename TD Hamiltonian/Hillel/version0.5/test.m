function test
global L
global x
global T
global w0
global matvecs
matvecs=0;
n=128;
[x,u0,L]=initial(n);
m1=8;
m2=8;
steps=50;
T=1;
dt=T/steps;
t=0;
w0=1;
tic
for k=1:steps
    [u,u0]=tdrhs_chb(u0,'matvec',t,dt,m1,m2);
    t=t+dt;
end
toc
plot(abs(u))
load fu1 u1
norm(u1-u)
matvecs

