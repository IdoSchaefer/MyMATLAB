function test
global A
n=50;
w=ones(n,1);
A=spdiags([w -2*w w],-1:1,n,n);
v=rand(n,1);
m=10;
[v,h]=arnoldi(v,'matvec',m);
hh=h(1:m,1:m);
eig(hh)

m