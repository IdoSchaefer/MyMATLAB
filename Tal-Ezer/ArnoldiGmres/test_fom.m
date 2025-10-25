function test_fom
global A
n=100;
w=ones(n,1);
A=spdiags([-w 2*w w],-1:1,n,n);
v=rand(n,1);
m=10;
b=ones(n,1);
tol=1.e-6;
cycles=150;
[x,r,res,matvecs,flag]=fom(b,'matvec',m,tol,cycles);
[res matvecs flag]
