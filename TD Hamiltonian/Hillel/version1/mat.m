function A=mat(matvec,n)
u=zeros(n,1);
for i=1:n
    u(i)=1;
    A(:,i)=feval(matvec,u);
    u(i)=0;
end