function U = test1
global thalf
%global mu
%global E0
global T
global H
global t
%mu=1;
%T=9000;
T=10;
%E0=8*pi/T;
%tol=1.e-7;
tol=1.e-13;
maxit=2000;
%u0=exact(0);
u0 = [1;0];
mcheb=14;
steps=100;
dt=T/steps;
er0 = zeros(1, steps);
U = zeros(2, steps+1);
U(:, 1) = u0;
%H=zeros(2);
H = diag([1 2])*(-1i);
iu=sqrt(-1);
t=0;
for k=1:steps
    thalf=t+dt/2;
%     H(1,2)=-iu*mu/2*E0*S(thalf);
%     H(2,1)=-iu*mu/2*E0*S(thalf);
    H(1, 2) = -1i*exp(1i*thalf);
    H(2, 1) = -1i*exp(-1i*thalf);
    [u,tvec,iter,flag]=chb_tdrhs_exact1(u0,H,dt,mcheb,tol,maxit);
    t=t+dt;
    u0=u(:,mcheb);
%     uex=exact(t);
%     er0(k,1)=norm(u0-uex)/norm(uex);
%     er(k,1)=abs(abs(uex(1))^2-abs(u0(1))^2);
%     [t er0(k,1)];
    P = u0.*conj(u0);
    er0(k) = (P(2) - sin(t)^2);
    U(:, k+1) = u0;
end
%[max(er0) max(er)]
max(abs(er0))
t
