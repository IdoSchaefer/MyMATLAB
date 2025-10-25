function test
global A
global mvecs
% A=problem1(n);
n=1000;
A=problem1(n);
[n,n]=size(A);
k=15;tol=1.e-4;
m=40;cycles=100000;tol=1.e-4;
eta=1;
ireal=1;
mvecs=0;
r=rand(n,1);
[v,eg,er,flag]=evarnoldi(r,'matvec',eta,n,k,m,tol,cycles,ireal);
toc
[mvecs flag]
for j=1:k
    res=matvec(v(:,j))-eg(j)*v(:,j);
    er0(j,1)=norm(res);
end
clear all

