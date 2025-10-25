function a=dvd_vec(x,f)
n=length(x);
a(:,1)=f(:,1);
w(:,1)=f(:,1);
for k=2:n
   w(:,k)=f(:,k);
   for j=k-1:-1:1
       w(:,j)=(w(:,j+1)-w(:,j))/(x(k)-x(j));
   end
   a(:,k)=w(:,1);
end
