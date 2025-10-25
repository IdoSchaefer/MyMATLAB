function R=maxdomain(matvec,n)
[v,h]=arnoldi(rand(n,1),matvec,10);
e=eig(h(1:10,1:10));
R=max(abs(e));