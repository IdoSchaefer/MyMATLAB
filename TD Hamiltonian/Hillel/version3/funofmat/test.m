% All the comments are mine (Ido).
function test
global A
global nvec
nvec=0;
% n=40;
% A=zeros(n);
% dx=6/n;
% d=10;
% for i=2:n-1
%     A(i,i-1)=1/dx^2-d/dx; 
%     A(i,i)=-2/dx^2;
%     A(i,i+1)=1/dx^2+d/dx;
% end
%  A(1,1)=-2/dx^2;A(1,2)=1/dx^2+d/dx;A(n,n-1)=1/dx^2-d/dx;A(n,n)=-2/dx^2;
 n=100;
% This function creates a matrix (A in this case) from the function
% (matvec1 in this case) that defines its operation on a vector.
A=mat('matvec1',n);
nvec=0;
u0=rand(n,1);
save fu0 u0
load fu0 u0
m=10;tol=1.e-6;cycles=60;
T=1;
k=1;
dt=T/k;
ireal=1;
u=u0;
for j=1:k
    u=newton(u,dt,'matvec','fun',m,tol,cycles,ireal);
end
 ex=expm(T*A)*u0;
% ex=funm(A,'fun')*u0;
er=norm(ex-u)/(1+norm(ex));
er
nvec
clear all



