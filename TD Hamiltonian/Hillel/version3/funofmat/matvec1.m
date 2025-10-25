function u=matvec1(v)
global nvec
nvec=nvec+1;
n2=length(v);
n=n2/2;
u(1:n,1)=v(n+1:n2,1);
for j=2:n
    u(n+j,1)=(v(j+1,1)-2*v(j,1)+v(j-1,1))*n^2;
end
u(n+1,1)=(v(2,1)-2*v(1,1))*n^2;
u(n2,1)=(-2*v(n,1)+v(n-1,1))*n^2;
