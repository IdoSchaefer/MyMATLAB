function [v,h]=arnoldi(v,matvec,m)
for j=1:m
   v(:,j+1)=feval(matvec,v(:,j));
   for i=1:j
      h(i,j)=v(:,i)'*v(:,j+1);
      v(:,j+1)=v(:,j+1)-h(i,j)*v(:,i);
   end
   h(j+1,j)=norm(v(:,j+1));
   if h(j+1,j) < 1.e-14
     h(j+1,j)=1.e-14;
   end
   v(:,j+1)=v(:,j+1)/h(j+1,j);
end
