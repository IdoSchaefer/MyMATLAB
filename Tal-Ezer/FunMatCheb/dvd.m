function [a,w]=dvd(w,x,f,ro)
n=length(x);
[m,m0]=size(w);
for k=m+1:n
   w(k,:)=f(k-m,:);
   for j=k-1:-1:1
       w(j,:)=ro*(w(j+1,:)-w(j,:))/(x(k)-x(j));
   end
   a(k-m,:)=w(1,:);
end
  
