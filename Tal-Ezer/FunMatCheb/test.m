function test
% All comments are mine.
global A
global nvec
global alpha beta
global x g gtag
% These are unused throughout the algorithm:
alpha=.9;beta=.9;
nvec=0;
n=16;
eps=1.e-3;
% This is an estimation of the p parameter (Eq. (3.13) in the draft of the
% new paper):
p=2*atan(eps^(1/n));
% x is the modified Chebyshev grid. g is dy/dx. gtag is (d^2)y/(dx)^2.
% g and gtag are required for the computation of the 2'nd derivative.
[x,g,gtag]=nonpolpoints(p,n);
% v is the initial state - the ground state with a unity momentum kick:
v=initial(n);
m=8;tol=1.e-8;restarts=2*600;
a=[];w=[];z=[];
% tpoints=[.2:.2:.8];
% One time point to be computed, t=1:
tpoints=1;
ireal=1;
iu=sqrt(-1);
parm(1)=0;
parm(2)=2;
parm(3)=2;
parm(4)=iu;
parm(5)=-iu*300;
% 
% 
%  A=mat('matvec2',n);
% % save fA A
% % load fA A
%  sort(eig(A))
tic
u=funofmat(v,tpoints,'matvec','fun',m,tol,restarts,parm);
toc
k=length(tpoints);
nvec
% Error check:
A=mat('matvec',n);
for j=1:k
    ex=expm(tpoints(j)*A)*v;
    er(j,1)=norm(ex-u(:,j))/(norm(ex)+1);
end
er
clear all





