function test
global thalf
global mu
global E0
global T
global H
global t
mu=1;
T=9000;
E0=8*pi/T;
tol=1.e-6;maxit=2000;
u0=exact(0);
mcheb=14;
steps=2;
dt=T/steps;
H=zeros(2);
iu=sqrt(-1);
t=0;
for k=1:steps
    thalf=t+dt/2;
    H(1,2)=-iu*mu/2*E0*S(thalf);
    H(2,1)=-iu*mu/2*E0*S(thalf);
    [u,tvec,iter,flag]=chb_tdrhs_exact(u0,H,dt,mcheb,tol,maxit);
    t=t+dt;
    u0=u(:,mcheb);
    uex=exact(t);
    er0(k,1)=norm(u0-uex)/norm(uex);
    er(k,1)=abs(abs(uex(1))^2-abs(u0(1))^2);
    [t er0(k,1)]
end
[max(er0) max(er)]
t
