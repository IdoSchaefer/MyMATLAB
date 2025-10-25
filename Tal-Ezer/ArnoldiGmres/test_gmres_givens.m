function test_gmres_givens
n=100;
m=10;
b=[1:n]';
tol=1.e-8;
cycles=150;
[x,res,matvecs,flag]=gmres_givens(b,'matvec',m,tol,cycles);
r=b-matvec(x);
exres=norm(r)/norm(b);
[exres res matvecs flag]
m

