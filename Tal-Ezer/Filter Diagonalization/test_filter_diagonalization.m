function test_filter_diagonalization
global A
n=500;
w=ones(n,1);
A=spdiags([w -2*w w],-1:1,n,n);
u=eig(full(A));
[uu,p]=sort(-u);
u=u(p);
tol=1.e-6;m=600;
R=4;
t=3000;
s=10;
v=chbeig(R,t,'matvec',n,m,s);
[v,r]=myqr(v);
h=(v'*(A*v))/(v'*v);
eg=eig(h);
[ee,p]=sort(-eg);
eg=eg(p);
for j=1:s
    ex(j,1)=u(j);
    er(j,1)=abs(ex(j)-eg(j));
end
[ex eg er]
  
