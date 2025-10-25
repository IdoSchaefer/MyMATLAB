function [x,r,res,matvecs,flag]=fom(b,matvec,m,tol,cycles)
flag=1;
normb=norm(b);
n=length(b);
x=zeros(n,1);
r=b;
matvecs=0;
for cycle=1:cycles
   normr=norm(r);
   v(:,1)=r/normr;
   for j=1:m
       v(:,j+1)=feval(matvec,v(:,j));
       matvecs=matvecs+1;
       for i=1:j
           h(i,j)=v(:,i)'*v(:,j+1);
           v(:,j+1)=v(:,j+1)-h(i,j)*v(:,i);
       end
       h(j+1,j)=norm(v(:,j+1));
       v(:,j+1)=v(:,j+1)/h(j+1,j);
   end
   e1=ones(m,1);
   z=normr*e1;
   y=h(1:m,1:m)\z;
   x=x+v(:,1:m)*y;
   r=b-feval(matvec,x);
%    error=abs(norm(r)*QT(m+1,1))/normb;
   res=norm(r)/normb
   if res < tol
       flag=0;
       return
    end
end