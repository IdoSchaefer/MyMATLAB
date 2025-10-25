function [x,res,matvecs,flag]=gmres_givens(b,matvec,m,tol,cycles)
flag=1;
normb=norm(b);
n=length(b);
x=zeros(n,1);
r=b;
matvecs=0;
for cycle=1:cycles
   rs(1,1)=1;
   c(1,1)=0;
   s(1,1)=0;
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
       rmat(1:j+1,j)=h(1:j+1,j);
       [rmat(:,j),rs,c,s]=givens(rmat(:,j),j,rs,c,s);
       res=abs(rs(j+1))*normr/normb;
        if res < tol
            yy=normr*trian(rmat,rs,j);
            x=x+v(:,1:j)*yy;
            r=b-v(:,1:j+1)*(h(1:j+1,1:j)*yy);
            flag=0;
            return
       end
   end
   yy=normr*trian(rmat,rs,j);
   x=x+v(:,1:m)*yy;
   r=b-feval(matvec,x);
end
               