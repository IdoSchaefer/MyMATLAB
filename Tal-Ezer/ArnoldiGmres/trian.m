function x=trian(a,b,n)
x(n,1)=b(n,1)/a(n,n);
for i=n-1:-1:1
   x(i,1)=(b(i,1)-a(i,i+1:n)*x(i+1:n,1))/a(i,i);
end