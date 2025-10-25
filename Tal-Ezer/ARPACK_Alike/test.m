function test
global A
n=100;
w=ones(n,1);
A=spdiags([ w -2*w w ],-1:1,n,n);
m=20;cycles=100;k=4;tol=1.e-6;
[v,eg,error,matvecs,flag]=evarnoldi('matvec',n,k,m,tol,cycles);
[matvecs flag]
for j=1:k
    res=matvec(v(:,j))-eg(j)*v(:,j);
    er=norm(res);
    [eg(j) er]
end

