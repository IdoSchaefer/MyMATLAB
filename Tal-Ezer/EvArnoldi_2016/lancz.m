function [v,h]=lancz(v,h,matvec,n,m,j)
   if j>1
       [s1,s2]=size(h);
       sh=h(1:s2,1:s2);
       sv=v(:,1:s2);
       r=h(s1,s2)*v(:,s2+1);
   else
       r=v;
       sv=[];
       sh=[];
   end
   [m1,m2]=size(sh);
   for i=1:m1
       for j=1:m2
           if abs(i-j)>1
               sh(i,j)=0;
           end
       end
   end
   [v,h,r,ierr]=mlanpro(matvec,n,m,r,sv,sh,1);
   h(m+1,m)=norm(r);
   v(:,m+1)=r/h(m+1,m);
   h=full(h);