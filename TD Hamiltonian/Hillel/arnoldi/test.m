function test
global L
global x
global T
global w0
global mmcheb
global matvecs
global thalf
mcheb1=14;
mcheb2=4;
matvecs=0;
n=128;
[x,u0,L]=initial(n);
mmcheb=mcheb1;
tol=1.e-10;
steps=100;
T=10;
dt=T/steps;
t=0;
maxit=100;
w0=1;
tic
for k=1:steps
    thalf=t+dt/2;
    [u,flag]=tdrhs(u0,'matvec',t,dt,mcheb1,mcheb2,tol,maxit);
    t=t+dt;
    u0=u(:,mcheb1);
end
toc
plot(abs(u0))
load fu1 u1
e=norm(u1-u0);
[mcheb1 mcheb2 matvecs e]

