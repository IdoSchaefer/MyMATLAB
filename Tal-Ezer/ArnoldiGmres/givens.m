function [w,rs,c,s]=givens(w,n,rs,c,s);
ep=1.d-16;
if n > 1
   for i=2:n
      cs=[c(i-1) s(i-1);-s(i-1) c(i-1)];
      w(i-1:i,1)=cs*w(i-1:i,1);
   end
end
gama=norm(w(n:n+1));
if gama==0
   gama=ep;
end
c(n,1)=w(n)'/gama;
s(n,1)=w(n+1)'/gama;
w(n,1)=gama;
w(n+1,1)=0;
rs(n+1,1)=-s(n,1)*rs(n,1);
rs(n,1)=c(n,1)*rs(n,1);
