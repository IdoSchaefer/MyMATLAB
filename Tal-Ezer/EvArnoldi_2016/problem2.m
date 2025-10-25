function A=problem2(n)
nn=2*n+1;
w=1;
for i=2:nn
     j=i-n-1;
    q(i-1,1)=-w*sqrt(n*(n+1)-j*(j-1));
end
nn=2*n+1;
q(nn)=0;
for i=1:nn
    q1(i,1)=(2/n)*i^2/2;
end
A=spdiags([q q1],-1:0,nn,nn);
A=A+A';
