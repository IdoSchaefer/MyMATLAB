function A=mat(matvec,n)
u=zeros(n,1);
for j=1:n
    u(j)=1;
    A(:,j)=feval(matvec,u);
    u(j)=0;
end